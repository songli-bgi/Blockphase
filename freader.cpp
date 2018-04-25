#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <zlib.h>
#include "freader.h"

using std::string;

#define FREADER_BUF_SIZE 64 * 1024 * 1024

freader::freader() : fname_(), fp_(NULL), buffer_(NULL), buffer_size_(FREADER_BUF_SIZE),
begin_(0), end_(0), line_(NULL), m_line_(0), is_open_(0), is_eof_(0) {}

bool freader::open(const char *fname)
{
	fname_ = fname;
	if (fname_.empty()) {
		fprintf(stderr, "%s line %d: fname cannot be empty\n", __FILE__, __LINE__);
		return 0;
	}

	// check file format, to be fix

	// open plain text and gzip file
	fp_ = (fname_.compare("-")) ? gzopen(fname_.c_str(), "rb") : gzdopen(fileno(stdin), "r");
	if (fp_ == NULL) {
		fprintf(stderr, "%s line %d: fail to open %s\n", __FILE__, __LINE__, fname_.c_str());
		return 0;
	}

	buffer_ = (char *) calloc(buffer_size_, sizeof(char));
	if (buffer_ == NULL) {
		fprintf(stderr, "%s line %d: fail to allocate memory\n", __FILE__, __LINE__);
		return 0;
	}
	begin_ = end_ = 0;

	m_line_ = 256;
	line_ = (char *) calloc(m_line_, sizeof(char));
	if (line_ == NULL) {
		fprintf(stderr, "%s line %d: fail to allocate memory\n", __FILE__, __LINE__);
		return 0;
	}

	is_open_ = 1;
	is_eof_ = 0;
	return 1;
}

void freader::close()
{
	if (fp_) gzclose(fp_);
	if (buffer_) free(buffer_);
	if (line_) free(line_);
}

// return -1 at error or EOF
int freader::getchar()
{
	if (begin_ < end_) return buffer_[begin_++];
	else { // has to refresh buffer
		if (is_eof_) return -1; // end of file

		// refresh buffer
		begin_ = 0;
		end_ = gzread(fp_, buffer_, buffer_size_);
		if (end_ < 0) {
			fprintf(stderr, "%s line %d: gzread fail\n", __FILE__, __LINE__);
			return -1;
		}
		if (end_ < buffer_size_) is_eof_ = 1;

		if (end_ > 0) {
			return buffer_[begin_++];
		} else return -1; // very special case, reach end of file
	}
	return -1;
}

// return NULL at EOF
const char *freader::getline(char line_delimiter)
{
	if (eof()) return NULL;

	bool found = 0;
	int c = 0;
	int l = 0;
	while ((c = getchar()) >= 0) {
		if (c != line_delimiter) {
			if (l + 1 >= m_line_) { // has to expand memory
				m_line_ += 256; // increase 256 bytes
				line_ = (char *) realloc(line_, m_line_);
			}
			line_[l++] = c;
		} else {
			found = 1; break;
		}
	}
	line_[l] = '\0';
	return line_;
}

// return 0 at EOF
int freader::read(void *buf, int n)
{
	if (eof()) return 0;
	char *p = (char *) buf;
	int c = 0;
	int i = 0;
	for (i = 0; i < n; i++) {
		c = getchar();
		if (c < 0) break;
		if (p != NULL) p[i] = c;
	}
	return i;
}

bool freader::eof()
{
	if (begin_ >= end_ && is_eof_) return 1;
	return 0;
}

void freader::set_buffer_size(int size)
{
	if (is_open_) {
		fprintf(stderr, "%s line %d: cannot call set_buffer_size if a file is opened\n", __FILE__, __LINE__);
	} else buffer_size_ = size > 0 ? size : FREADER_BUF_SIZE;
}

