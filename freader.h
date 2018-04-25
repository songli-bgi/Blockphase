/**
 * @file
 * @brief  A generic file reader.
 * @date   24-Oct-2012, First realse.
 * @date   01-Nov-2012, Add drop feature to read().
 *
 * Supported file formats: plain text and gzip.
 *
 * More file formats will be supported in the futher.
 */

#include <string>
#include <stdint.h>
#include <zlib.h>

using std::string;

#ifndef __GZREADER_H__
#define __GZREADER_H__

/**
 * Use this class to read input files in multiple formats.
 * A intermediate buffer layer is applied to reduce IO burden 
 * Currently, plain text and gzip files are supported.
 * More file formats will be supported in the futher.
*/
class freader
{
public:
	/// Create a file reader.
	freader();

	/**
	 * Open a input file.
	 * @param fname Name of a input file.
	 * @return \c true if the input file is successfully opened, \c false if fail.
	 * Supported file formats: plain text and gzip.
	 * To read from STDIN, passing "-" to fname.
	 */
	bool open(const char *fname);

	/// Close a opened file.
	void close();

	/**
	 * Get a character from the input file.
	 * @return A character if success, -1 at end-of-file.
	 */
	int getchar();

	/**
	 * Get a line from the input file.
	 * @param line_delimiter Line delimiter, default is '\\n'.
	 * @return A new line as a c-string. Line delimiter is not included.
	 * @return NULL at EOF.
	 * @note If the input file is compressed, it will be auto-decompressed.
	 */
	const char *getline(char line_delimiter = '\n');

	/**
	 * Read n decompressed bytes from the input file.
	 * @param buf The location to store the n bytes reading from the input file.
	 * @param n   The number of bytes to read out.
	 * @return The number of uncompressed bytes actually read, 0 and a short number indicates end of file, 0 at EOF.
	 * @note If \c buf is NULL, the bytes read out with be dropped.
	 */
	int read(void *buf, int n);

	/**
	 * Test end-of-file.
	 * @return \c true if reach end-of-file, \c false if not.
	 */
	bool eof();

	/**
	 * Set buffer size.
	 * @param size New buffer size, default is 64 * 1024.
	 * This method has to be called before openning a file.
	 */
	void set_buffer_size(int size);

private:
	// prevent copy
	freader(const freader &);

	// prevent assignment
	freader & operator= (const freader &);

	string fname_;
	gzFile fp_;
	char *buffer_;
	int32_t buffer_size_;  // default is 64 * 1024
	int32_t begin_;        // valid begin position in the buffer
	int32_t end_;          // one off the end position of the buffer
	char *line_;           // keep the current line
	int32_t m_line_;       // memory usage of _line
	bool is_open_;    // opened or not
	bool is_eof_;     // reach EOF or not
};

#endif

/**
\example
\code
#include <stdio.h>
#include "freader.h"

int main (int argc, char **argv)
{
	// create a freader
	freader inf;
	int ret;

	// open
	if (argc >= 2) {
		ret = inf.open(argv[1]);
		if (ret == false) {
			fprintf(stderr, "fail to open %s\n", argv[1]);
			return 0;
		}
	} else { // read from stdin
		ret = inf.open("-");
		if (ret == false) {
			fprintf(stderr, "fail to open stdin\n");
			return 0;
		}
	}

	// getchar
	int c;
	while ((c = inf.getchar()) >= 0) putchar(c);

	// getline
	const char *s;
	while ((s = inf.getline()) != NULL)
		printf("%s\n",s);

	inf.close();
	return 0;
}
\endcode
*/

