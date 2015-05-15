#include <xstream/config.h>
#include <xstream/digest.h>
#include <streambuf>
#include <algorithm>

#include "debug.h"


namespace xstream{
namespace digest{

	static const int eof = std::streambuf::traits_type::eof();

	stream::stream(size_t sz)
		:buf(sz),length(0)
	{
		LOG("digest::stream");
		setp (buf.buf, buf.buf + buf.size);
	}

	stream::~stream()
	{
		LOG("digest::~stream");
	}

	int stream::sync()
	{
		LOG("digest::stream::sync");

		calculate_digest();
		setp (buf.buf, buf.buf + buf.size);
		return 1;
	}

	int stream::overflow(const int c)
	{
		LOG("digest::stream::overflow "<<c);

		if ( eof == c ) {
			LOG ("\tEOF");
			//XXX what should I do?
			return eof;
		}else{
			sync();
			*(pptr()) = static_cast < char >(c);
			pbump (1);
		}

		pbump (1);
		length+=buf.size;
		LOG("\tlength = "<<length);
		return 1;
	}

	std::streamsize stream::xsputn(char *b, std::streamsize n)
	{
		LOG("digest::stream::xsputn("<<b<<", "<<n<<")");

		sync();
		setp (b, b + n);
		length+=n;
		LOG("\tlength = "<<length);
		sync();
		
		return n;
	}

	block_stream::block_stream(const size_t chunk, const unsigned int blocks)
		:stream(chunk*blocks),chunk_size(chunk)
	{
		LOG("digest::block_stream ("<<chunk<<", "<<blocks<<")");
	}

	int block_stream::sync()
	{
		LOG("digest::block_stream::sync");
		
		//number of bytes in the buffer
		size_t t = taken();
		char* beg = buf.buf;

		while( t>=chunk_size ){
			LOG("\tchunk");
			setp(beg, beg+chunk_size);
			calculate_digest();
			beg+=chunk_size;
			t-=chunk_size;
			length+=chunk_size;
			LOG("\tlength = "<<length);
		}

		//reset buffer
		
		if(0==t){
			LOG("\tnothing remains, reset buffer");
			setp (buf.buf, buf.buf + buf.size);
		}else if(beg != buf.buf){
			LOG("\t"<<t<<" bytes remain, copying");
			//t bytes at tail, need to be placed at the beggining and buffer set acordingly
			std::copy(beg, beg+t, buf.buf);
			setp (buf.buf, buf.buf + buf.size);
			pbump(t);
		}

		return 1;
	}

	int block_stream::overflow(const int c)
	{
		LOG("digest::block_stream::overflow "<<c);

		if ( eof == c ) {
			LOG ("\tEOF");
			//XXX what should I do?
			return eof;
		}else{
			sync();
			*(pptr()) = static_cast < char >(c);
			pbump (1);
		}

		return 1;
	}


	std::streamsize block_stream::xsputn(char *b, std::streamsize n)
	{
		LOG("digest::block_stream::xsputn "<<n);
		sync();
		const std::streamsize t = taken();
		const std::streamsize r = chunk_size - t;
		if(t>0){
			LOG("\t"<<t<<" bytes previously in buffer, copying first "<<r<<" bytes");
			std::copy(b, b+r, buf.buf+t);
			sync();
		}
		setp(b+r, b+n);
		sync();

		return n;
	}

}//namespace digest
}//namespace xstream
