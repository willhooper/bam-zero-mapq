use std::ffi::CString;
use std::io;
use std::path::Path;
use std::ptr::NonNull;

use rust_htslib::bam::{HeaderView, Record};
use rust_htslib::htslib;

use crate::contextual;

// rust-htslib 0.44.1 hides the Writer handle and the ThreadPool handle. Use its
// HTSlib bindings for output so sam_idx_init/sam_idx_save can build the index
// during writing while input and output still share one compression pool.
pub struct ThreadPool {
    pool: NonNull<htslib::hts_tpool>,
}

impl ThreadPool {
    pub fn new(threads: u32) -> io::Result<Self> {
        let threads = i32::try_from(threads)
            .ok()
            .filter(|threads| *threads > 0)
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "thread count is out of range")
            })?;
        // SAFETY: HTSlib accepts a positive thread count and returns an owned pool.
        let pool = NonNull::new(unsafe { htslib::hts_tpool_init(threads) }).ok_or_else(|| {
            contextual(
                "HTSlib thread pool allocation failed",
                io::Error::last_os_error(),
            )
        })?;
        Ok(Self { pool })
    }

    /// # Safety
    /// The file must be open and must be closed before this pool is dropped.
    pub unsafe fn attach(&self, file: *mut htslib::htsFile) -> io::Result<()> {
        let mut config = htslib::htsThreadPool {
            pool: self.pool.as_ptr(),
            qsize: 0, // Let HTSlib choose its default queue size.
        };
        // HTSlib copies the configuration; it borrows the pool itself.
        if unsafe { htslib::hts_set_thread_pool(file, &mut config) } != 0 {
            return Err(contextual(
                "HTSlib thread pool attachment failed",
                io::Error::last_os_error(),
            ));
        }
        Ok(())
    }
}

impl Drop for ThreadPool {
    fn drop(&mut self) {
        // SAFETY: all attached files must have been closed before destroying the pool.
        unsafe { htslib::hts_tpool_destroy(self.pool.as_ptr()) };
    }
}

pub struct IndexedBamWriter<'pool> {
    file: Option<NonNull<htslib::htsFile>>,
    header: HeaderView,
    // sam_idx_init borrows this string until sam_idx_save; keep it through close.
    index_path: CString,
    _pool: &'pool ThreadPool,
}

impl<'pool> IndexedBamWriter<'pool> {
    pub fn new(
        path: &Path,
        header: HeaderView,
        compression: u32,
        pool: &'pool ThreadPool,
    ) -> io::Result<Self> {
        let mut index_path = path.as_os_str().to_os_string();
        index_path.push(".bai");
        let index_path = CString::new(index_path.as_encoded_bytes())
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error))?;
        let path = CString::new(path.as_os_str().as_encoded_bytes())
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidInput, error))?;
        let mode = CString::new(format!("wb{compression}")).expect("mode contains no NUL");

        // SAFETY: both strings are NUL-terminated; hts_open returns an owned handle.
        let file = NonNull::new(unsafe { htslib::hts_open(path.as_ptr(), mode.as_ptr()) })
            .ok_or_else(|| contextual("cannot open BAM", io::Error::last_os_error()))?;
        let mut writer = Self {
            file: Some(file),
            header,
            index_path,
            _pool: pool,
        };

        // SAFETY: writer owns the open handle and header and borrows the pool for
        // its full lifetime. Drop closes the handle on initialization errors.
        unsafe {
            pool.attach(file.as_ptr())?;
            if htslib::sam_hdr_write(file.as_ptr(), writer.header.inner_mut()) < 0 {
                return Err(contextual(
                    "cannot write BAM header",
                    io::Error::last_os_error(),
                ));
            }
            // Must run after the header and before the first record. min_shift=0
            // selects BAI. sam_write1 then feeds index entries to HTSlib's BGZF
            // writer, which resolves offsets as compressed blocks are emitted.
            if htslib::sam_idx_init(
                file.as_ptr(),
                writer.header.inner_mut(),
                0,
                writer.index_path.as_ptr(),
            ) < 0
            {
                return Err(contextual(
                    "cannot initialize BAI index",
                    io::Error::last_os_error(),
                ));
            }
        }
        Ok(writer)
    }

    pub fn write(&mut self, record: &Record) -> io::Result<()> {
        let file = self.file.expect("writer is open");
        // SAFETY: the file, header, and record remain valid throughout the call.
        // sam_write1 is required here: bam_write1 would bypass index updates.
        if unsafe { htslib::sam_write1(file.as_ptr(), self.header.inner(), record.inner()) } < 0 {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "BAM write/index update failed; input must be coordinate-sorted and within BAI limits",
            ));
        }
        Ok(())
    }

    pub fn finish(mut self) -> io::Result<()> {
        let file = self.file.expect("writer is open");
        // SAFETY: all records have been submitted and the file, pool, header, and
        // index filename are alive. sam_idx_save drains BGZF workers and saves
        // the accumulated index; it does not reopen or reread the BAM.
        if unsafe { htslib::sam_idx_save(file.as_ptr()) } < 0 {
            return Err(contextual(
                "cannot finish BAI index (check sort order, BAI coordinate limits, and output storage)",
                io::Error::last_os_error(),
            ));
        }

        // hts_close consumes the handle even on failure. Take it first so Drop
        // cannot close it twice, and report delayed compression/write errors.
        let file = self.file.take().expect("writer is open");
        if unsafe { htslib::hts_close(file.as_ptr()) } < 0 {
            return Err(contextual("cannot close BAM", io::Error::last_os_error()));
        }
        Ok(())
    }
}

impl Drop for IndexedBamWriter<'_> {
    fn drop(&mut self) {
        if let Some(file) = self.file.take() {
            // SAFETY: this is the sole owner; header, filename, and pool are alive.
            // Error paths close the BAM without saving an incomplete index.
            unsafe { htslib::hts_close(file.as_ptr()) };
        }
    }
}
