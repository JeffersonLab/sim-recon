#include <xstream/config.h>

#if HAVE_LIBBZ2

#include <algorithm>
#include <cassert>

#include <xstream/bz.h>
#include <xstream/except/bz.h>

#include <bzlib.h>


#include "debug.h"

namespace xstream {
namespace bz {

    static const int eof = std::streambuf::traits_type::eof();

    struct pimpl: public bz_stream {};

    static inline int flush_macro(const flush_kind f) {
		switch (f) {
			case no_sync:
				return BZ_RUN;
			case full_sync:
				return BZ_FLUSH;
			case finish_sync:
				return BZ_FINISH;
			default:
				//should probably throw
				return BZ_RUN;
		}
	}

	common::common(std::streambuf * sb)
	: xstream::common_buffer (sb), z_strm (0) {
		LOG("bz::common");	

		z_strm = new pimpl;

		//initialize zlib structure
		z_strm->bzalloc = NULL;
		z_strm->bzfree = NULL;
		z_strm->opaque = NULL;
		//buffers
		z_strm->avail_out = out.size;
		z_strm->next_out = out.buf;

		z_strm->avail_in = 0;
		z_strm->next_in = in.buf;

	}

	unsigned long int common::input_count() const {
		return z_strm->total_in_lo32;
	}

	unsigned long int common::output_count() const {
		return z_strm->total_out_lo32;
	}

	common::~common() {
		LOG("bz::~common");
		delete z_strm;
	}


	//default compression 9
	ostreambuf::ostreambuf (std::streambuf * sb):common(sb),level (9) {
		LOG ("bz::ostreambuf without compression level");
		init ();
	}

	ostreambuf::ostreambuf (std::streambuf * sb, const int l):
		common(sb),level (l) {
			LOG ("bz::ostreambuf with compression level " << l);
			init ();
		}

	const char* error_str(const int err){
		switch(err){
			case BZ_MEM_ERROR:
				return "out of memory";
			case BZ_CONFIG_ERROR:
				return "bzlib badly configured (bad sizes of int32 (4), int16 (2) or char (1), check and recompile)";
			case BZ_PARAM_ERROR:
				return "invalid parameter, possibly invalid compression level";
			case BZ_SEQUENCE_ERROR:
				return "bad sequence (this means xstream is buggy)";
			case BZ_DATA_ERROR:
				return "invalid or incomplete data (crc failed)";
			case BZ_DATA_ERROR_MAGIC:
				return "magic bytes not found in stream";
			case BZ_IO_ERROR:
				return "io error";
			case BZ_UNEXPECTED_EOF:
				return "premature end of data";
			case BZ_OUTBUFF_FULL:
				return "output buffer full";
		}
		
		return "unknown error";
	}

	void ostreambuf::raise_error(const int err){
		std::string what = error_str(err);

		LOG("bz::ostreambuf::raise_error ("<<err<<") = "<<what);

		if(what.size()>0){
			throw compress_error(this,what);
		}else{
			throw compress_error(this);
		}

	}


	void ostreambuf::init () {
		LOG ("bz::ostreambuf::init");
		int cret =::BZ2_bzCompressInit(
			z_strm,
			level, 
			0, //verbosity
			0  //workFactor (default value) controls when to switch to the fallback algorithm
		);

		if (BZ_OK != cret) {
			LOG ("bz::ostreambuf::init: error creating zstream " << cret);
			raise_error(cret);
		}
		//initialize streambuf interface functions
		setp (in.buf, in.buf + in.size);
	}

	ostreambuf::~ostreambuf () {
		LOG ("bz::ostreambuf::~ostreambuf");
		//fullsync (write remaining data)
		flush (finish_sync);

		//sync underlying streambuf
		_sb->pubsync();

		if (0 != z_strm) {
			//XXX should I throw an exception in case of error?
			//remember this is a destructor
			int cret = ::BZ2_bzCompressEnd(z_strm);
			if(BZ_OK != cret){
				LOG("\tERROR: BZ2_bzCompressEnd returned "<<cret);
			}
		}
	}

	int ostreambuf::sync () {
	LOG ("bz::ostreambuf::sync");
	int ret = flush (full_sync);
	  _sb->pubsync();
	  return ret;
	}


	int ostreambuf::overflow (const int c) {
		LOG ("bz::ostreambuf::overflow(" << c << ")\t available=" << (available ()) << "\tEOF=" << eof);
		if ( eof == c ) {
			LOG ("\tEOF");
			//XXX what should I do?
			return eof;
		} else {
			if (0 == available ()) {
				LOG ("\t have to flush :[" << in.buf << "]");
				flush (no_sync);
			}
			*pptr () = static_cast < char >(c);

			pbump (1);
		}
		return c;
	}

	std::streamsize ostreambuf::xsputn (char *buffer, std::streamsize n) {
		LOG ("bz::ostreambuf::xsputn(" << buffer << "," << n << ")");

		std::streamsize written = z_strm->avail_in;

		flush(no_sync);
		//nothing should be on input buffer
		assert( 0==z_strm->avail_in );

		//store original zlib stream state
		size_t real_s_in  = in.size;
		char* real_in = in.buf;

		//this is very tricky
		//the destructor cannot be called during this block
		//because that would cause a free to the given buffer
		try{
			//fake that the buffer is the new input buffer
			in.size = n;
			in.buf  = buffer;

			flush(no_sync);
		}
		catch(...){
			in.size = real_s_in;
			in.buf  = real_in;
			throw;
		}

		//restore zlib stream state

		in.size = real_s_in;
		in.buf  = real_in;

		written += z_strm->next_in - in.buf;

		
		return written;
	}

	int ostreambuf::flush (flush_kind f) {
		LOG ("bz::ostreambuf::flush(" << f << ")");
		std::streamsize in_s = taken ();
		LOG ("\tinput_size=" << in_s);

		//reset input
		z_strm->next_in = pbase();
		z_strm->avail_in = in_s;

		bool redo = false;

		do {
			int cret;
			redo=false;
			bool error=false;

			//reset output
			z_strm->next_out = out.buf;
			z_strm->avail_out = out.size;

			/*
			 LOG ("\tpre_deflate: [" << z_strm->next_in << "]\tavail_in=" << z_strm->avail_in <<
					"\tavail_out=" << z_strm->avail_out);
			*/

			cret =::BZ2_bzCompress (z_strm, flush_macro(f) );

			/*
			 LOG ("\tpost_deflate: [" << z_strm->next_in << "]\tavail_in=" << z_strm->avail_in <<
					"\tavail_out=" << z_strm->avail_out << "\tcret=" << cret);
			*/

			//error handling
			if ( finish_sync==f ) {
				if (BZ_STREAM_END == cret) {
					redo=false;
				//XXX manual mentions BZ_FINISHING but this macro isn't defined anywhere in the source code of the library
				// so I use BZ_FINISH_OK but I'm not sure
				}else if(BZ_FINISH_OK == cret) {
					redo=true;
				}else{
					//serious error, throw exception
					LOG ("\terror in finish:" << cret);
					error=true;
				}
			} else if( full_sync==f ){
				if(BZ_FLUSH_OK == cret){
					LOG("\tanother go at sync");
					redo=true;
				}else if(BZ_RUN_OK==cret){
					LOG("\tsync ok");
					redo=false;
				}else{
					LOG("\terror in sync: "<<cret);
					error=true;
				}
			} else if(no_sync==f){
				if (BZ_RUN_OK != cret) {
					LOG ("\terror compressing " << cret);
					error=true;
				}
			}else{
				LOG("\tERROR: unkown flush mode "<<flush_macro(f));
				throw general_error();
			}

			if(error){
				raise_error(cret);
			}

			LOG ("\twritting " << (out.size - z_strm->avail_out) << " bytes");

			const std::streamsize count = out.size - z_strm->avail_out;
			const std::streamsize wrote = _sb->sputn (out.buf, count);
			if(count != wrote){
				LOG("\terror writting, only wrote "<<wrote<<" but asked for "<<count);
				raise_error(BZ_IO_ERROR);
			}

			if (0 == z_strm->avail_out){ // && 0 != z_strm->avail_in)
				LOG("\tavail_out=0 => redo");
				redo = true;
			}
		}
		while (redo);
		assert (0 == z_strm->avail_in);
		//reset buffer
		setp(in.buf, in.buf + in.size);
		return 0;

	}

	/////////////////////
	// istream follows //
	/////////////////////

	istreambuf::istreambuf (std::streambuf * sb):common(sb), end(false) {
		LOG ("bz::istreambuf");

		int cret =::BZ2_bzDecompressInit(z_strm,
			0, //verbosity
			0  //no small memory
		);

		if (BZ_OK != cret) {
			LOG ("\terror creating zstream " << cret);
			raise_error(cret);
		}
		//initialize streambuf interface functions
		//first call will call uflow and this will set the buffer accordingly
		//no buffering
		setg(out.buf,out.buf,out.buf);
	}

	void istreambuf::raise_error(const int err){
		std::string what = error_str(err);

		LOG("bz::istreambuf::raise_error ("<<err<<") = "<<what);

		if(what.size()>0){
			throw decompress_error(this,what);
		}else{
			throw decompress_error(this);
		}
	}

	int istreambuf::underflow() {
		LOG("z:istreambuf::underflow");

		if(end){
			LOG("\tend of stream (EOF)");
			//signal the stream has reached it's end
			return eof;
		}

		z_strm->avail_out = out.size;
		z_strm->next_out = out.buf;

		if(0 < z_strm->avail_in ){
			LOG("\tdata in queue, inflating");
			decompress();
		}

		while( !(end || 0==z_strm->avail_out) ){
			read_decompress();
		}
			
		//set streambuf pointers
		setg(out.buf, out.buf, z_strm->next_out);

		return 0;
	}
	
	//read to buffer in place (apart from data already buffered)
	std::streamsize istreambuf::xsgetn(char *buffer, std::streamsize n) {
		LOG("bz::istreambuf::xsgetn ("<<n<<")");

		if(end){
			LOG("\tend of stream (EOF)");
			//signal the stream has reached it's end
			return eof;
		}

		//copy input buffered data to buffer

		std::streamsize read = pbase() - pptr();
		std::copy(pbase(), pptr(), buffer);

		//store original zlib stream state
		size_t real_s_out = out.size;
		char*  real_out   = out.buf;

		try{
			out.buf  = buffer + read;
			out.size = n - read;;
			uflow();
		}catch(...){
			out.buf   = real_out;
			out.size = real_s_out;
			throw;
		}

		out.buf  = real_out;
		out.size = real_s_out;
		
		//next read will call uflow
		setg(out.buf,out.buf,out.buf);
		return (z_strm->next_out - buffer);
	}

	void istreambuf::read_decompress() {
		LOG("bz::istreambuf::read_inflate ");

		z_strm->next_in = in.buf;
		
		int read = _sb->sgetn(in.buf, in.size);
		LOG("\tread "<<read<<" bytes");

		if(0==read){
			LOG("\tpremature end of stream");
			raise_error(BZ_UNEXPECTED_EOF);
		}

		z_strm->avail_in=read;

		decompress();

	}

	void istreambuf::decompress() {
		LOG("bz::istreambuf::inflate ");

		int cret = ::BZ2_bzDecompress(z_strm);

		if(BZ_STREAM_END == cret){
			end=true;
		}else if(BZ_OK != cret){
			LOG("\terror inflating: "<<cret);
			raise_error(cret);
			//can try to salvage some more data with inflateSync (on some cases)
		}

	}

	/*XXX Should add a sync method */

	istreambuf::~istreambuf() {
		LOG("bz::~istreambuf");
		if (0 != z_strm) {
			int cret = ::BZ2_bzDecompressEnd(z_strm);
			if(BZ_OK != cret){
				LOG("\tERROR: BZ2_bzDecompressEnd returned "<<cret);
			}
		}
	}

}//namespace bz
}//namespace xstream

#endif //bzlib
