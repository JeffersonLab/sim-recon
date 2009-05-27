#include <xstream/config.h>

#if HAVE_LIBZ

#include <algorithm>
#include <string>

#include <xstream/z.h>
#include <xstream/except/z.h>
#include <stdexcept>

#include <zlib.h>

#include <cassert>

#include "debug.h"

namespace xstream {
namespace z {

    struct pimpl: public z_stream {};
	
    static const int eof = std::streambuf::traits_type::eof();

    static inline int flush_macro(const flush_kind f) {
		switch (f) {
			case no_sync:
				return Z_NO_FLUSH;
			case normal_sync:
				return Z_SYNC_FLUSH;
			case full_sync:
				return Z_FULL_FLUSH;
			case finish_sync:
				return Z_FINISH;
			case block_sync:
#ifndef Z_BLOCK
# define Z_BLOCK 5
#endif
				return Z_BLOCK;
			default:
				//should probably throw
				return Z_NO_FLUSH;
		}
	}

	const char* error_str(const int err){
		switch(err){
			case Z_MEM_ERROR:
				return "out of memory";
			case Z_VERSION_ERROR:
				return "zlib version mismatch";
			case Z_DATA_ERROR:
				return "invalid or incomplete data";
			case Z_STREAM_ERROR:
				return "stream error";
			case Z_NEED_DICT:
				return "need dictionary";
			case Z_STREAM_END:
				return "stream end";
			case Z_BUF_ERROR:
				return "buffer error";
		}

		return "unknown error";
	}
	

	common::common(std::streambuf * sb)
	: xstream::common_buffer (sb), z_strm (0) {
			LOG("z::common");	

			z_strm = new pimpl;

			//initialize zlib structure
			z_strm->zalloc = Z_NULL;
			z_strm->zfree = Z_NULL;
			z_strm->opaque = Z_NULL;
			//buffers
			z_strm->avail_out = out.size;
			z_strm->next_out = reinterpret_cast < Bytef* >(out.buf);

			z_strm->avail_in = 0;
			z_strm->next_in = reinterpret_cast < Bytef* >(in.buf);

		}

	void common::grow_out (const unsigned int factor) {

		const size_t taken = out.size - z_strm->avail_out;

		out.grow(factor);

		z_strm->next_out = reinterpret_cast < Bytef* >(out.buf + taken);
		z_strm->avail_out = out.size - taken;
	}

	unsigned long int common::input_count() const {
		return z_strm->total_in;
	}

	unsigned long int common::output_count() const {
		return z_strm->total_out;
	}

	unsigned long int common::checksum() const {
		return z_strm->adler;
	}

	common::~common() {
		LOG("z::~common");
		delete z_strm;
	}


	ostreambuf::ostreambuf (std::streambuf * sb):common(sb),level (Z_DEFAULT_COMPRESSION) {
		LOG ("z::ostreambuf without compression level");
		init ();
	}

	ostreambuf::ostreambuf (std::streambuf * sb, const int l):
		common(sb),level (l) {
			LOG ("z::ostreambuf with compression level " << l);
			init ();
		}

	void ostreambuf::raise_error(const int err) {
		std::string what = error_str(err);

		LOG("z::ostreambuf::raise_error ("<<err<<") = "<<what);

		if(what.size()>0){
			throw compress_error(this,what);
		}else{
			throw compress_error(this);
		}

	}


	void ostreambuf::init () {
		LOG ("z::ostreambuf::init");

		if(Z_DEFAULT_COMPRESSION==level || (level<=9 &&level>=1)){
			int cret =::deflateInit (z_strm, level);

			if (Z_OK != cret) {
				LOG ("z::ostreambuf::init: error creating zstream " << cret);
				//XXX exception ins constructor
				raise_error(cret);
			}
			//initialize streambuf interface functions
			setp (in.buf, in.buf + in.size);
		}else{
			throw std::domain_error("invalid compression level "+level);
		}

	}

	ostreambuf::~ostreambuf () {
		LOG ("z::ostreambuf::~ostreambuf");
		//sync (write remaining data)
		flush (finish_sync);

		//sync underlying streambuf
		_sb->pubsync();

		if (0 != z_strm) {
			//XXX should I throw an exception in case of error?
			//remember this is a destructor
			//I should definitly LOG something
			int cret = ::deflateEnd (z_strm);
			if(Z_OK!=cret){
				LOG("z::~ostreambuf error dealocating zstream");
			}
		}
	}

    int ostreambuf::sync () {
      LOG ("z::ostreambuf::sync");
      int ret = flush (normal_sync);
	  _sb->pubsync();
	  return ret;
    }


	int ostreambuf::overflow (const int c) {
		LOG ("z::ostreambuf::overflow(" << c << ")\t available=" << (available ()) << "\tEOF=" << eof);
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
		LOG ("z::ostreambuf::xsputn(" << buffer << "," << n << ")");

		//XXX should written be defined this way???
		std::streamsize written = z_strm->avail_in;

		flush(no_sync);
		//nothing should be on input buffer
		assert( 0==z_strm->avail_in );

		//store original zlib stream state
		char*  real_buf  = in.buf;
		size_t real_size = in.size;

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
			LOG("\t caught exception. restoring original values");
			in.buf  = real_buf;
			in.size = real_size;
			throw;
		}

		//restore zlib stream state
		in.buf  = real_buf;
		in.size = real_size;

		written += reinterpret_cast<char*>(z_strm->next_in) - in.buf;

		return written;
	}

	int ostreambuf::flush (const flush_kind f) {
		LOG ("z::ostreambuf::flush(" << f << ")");
		std::streamsize in_s = taken ();
		LOG ("\tinput_size=" << in_s);

		//reset input
		z_strm->next_in = reinterpret_cast < Bytef* >(pbase());
		z_strm->avail_in = in_s;

		bool redo = false;

		do {
			int cret;
			redo=false;

			//reset output
			z_strm->next_out = reinterpret_cast < Bytef* >(out.buf);
			z_strm->avail_out = out.size;

			do {
				cret =::deflate (z_strm, flush_macro (f));

				if(finish_sync==f && Z_OK==cret){
					grow_out();
					continue;
				}else{
					break;
				}
			} while (1);

			//error handling
			if (f == finish_sync) {
				if (Z_STREAM_END != cret) {
					//serious error, throw exception
					LOG ("\terror :" << cret);
				}
			} else if (Z_OK != cret) {
				LOG ("\terror deflating " << cret);
				//XXX throw exception here
				raise_error(cret);
			}

			LOG ("\twritting " << (out.size - z_strm->avail_out) << " bytes");
			//XXX need to check return value and wrap this
			_sb->sputn (out.buf, out.size - z_strm->avail_out);

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
		LOG ("z::istreambuf");

		int cret =::inflateInit (z_strm);

		if (Z_OK != cret) {
			LOG ("\terror creating zstream " << cret);
			//XXX throw exception here
			raise_error(cret);
		}
		//initialize streambuf interface functions
		//first call will call uflow and this will set the buffer accordingly
		//no buffering
		setg(out.buf,out.buf,out.buf);
	}

	void istreambuf::raise_error(const int err){
		std::string what = error_str(err);

		LOG("z::istreambuf::raise_error ("<<err<<") = "<<what);

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
		z_strm->next_out = reinterpret_cast < Bytef* >(out.buf);

		if(0 < z_strm->avail_in ){
			LOG("\tdata in queue, inflating");
			inflate();
		}

		while( !(end || 0==z_strm->avail_out) ){
			read_inflate();
		}
			
		//set streambuf pointers
		setg(out.buf, out.buf, reinterpret_cast <char*> (z_strm->next_out) );

		return 0;
	}
	
	//read to buffer in place (apart from data already buffered)
	std::streamsize istreambuf::xsgetn(char *buffer, std::streamsize n) {
		LOG("z::istreambuf::xsgetn ("<<n<<")");

		if(end){
			LOG("\tend of stream (EOF)");
			//signal the stream has reached it's end
			return eof;
		}

		//copy input buffered data to buffer

		std::streamsize read = pbase() - pptr();
		std::copy(pbase(), pptr(), buffer);

		//store original zlib stream state
		char* real_buf   = out.buf;
		size_t real_size = out.size;

		try{
			out.buf   = buffer + read;
			out.size  = n - read;;
			uflow();
		}catch(...){
			LOG("\t caught exception. restoring original values");
			out.buf  = real_buf;
			out.size = real_size;
			throw;
		}

		out.buf  = real_buf;
		out.size = real_size;
		
		//next read will call uflow
		setg(out.buf,out.buf,out.buf);
		return (reinterpret_cast<char*>(z_strm->next_out) - buffer);
	}

	void istreambuf::read_inflate( const flush_kind f) {
		LOG("z::istreambuf::read_inflate "<<f);

		z_strm->next_in = reinterpret_cast < Bytef* >(in.buf);
		
		int read = _sb->sgetn(in.buf,in.size);
		LOG("\tread "<<read<<" bytes");

		if(0==read){
			LOG("\tpremature end of stream");
			//XXX throw exception
			raise_error(Z_DATA_ERROR);
		}

		z_strm->avail_in=read;

		inflate(f);

	}

	void istreambuf::inflate(const flush_kind f) {
		LOG("z::istreambuf::inflate "<<f);

		int cret = ::inflate(z_strm, flush_macro(f));

		if(Z_STREAM_END == cret){
			end=true;
		}else if(Z_OK != cret){
			LOG("\terror inflating: "<<cret);
			//XXX throw exception
			raise_error(cret);
			//can try to salvage some more data with inflateSync (on some cases)
		}

	}

	/*XXX Should add a sync method */

	istreambuf::~istreambuf() {
		LOG("z::~istreambuf");
		if (0 != z_strm) {
			//XXX should I throw an exception in case of error?
			//remember this is a destructor
			::inflateEnd (z_strm);
		}
	}

}//namespace z
}//namespace xstream

#endif	//zlib
