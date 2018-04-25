/**
 * @file
 * @brief Parsing command line arguments.
 * @date  30-Oct-2012
 *
 * Features:
 * - Single letter options head with one hyphen '-' are parsed following POSIX convention.
 * - Long options head with two hyphens '--' are parsed following GNU conventions.
 * - options can mix with non-option arguments, although options always precede non-option arguments.
 * - Argument -- stops argument parsing. Arguments following -- are collected as non-option arguments.
 * - A single hypen - is interpreted as an ordinary non-option arguments.
 * - Ordinary non-option arguments are collected for the further use.
 * - Mixup of short and long options at command line is supported.
 * - Specifying a option more than once is not supported.
 *
 * In the writing of command line, short options should follow POSIX convention
 * and long options should follow the GNU extensions. E.g, suppose program foo
 * has two flags (-a -all, -v --verbose) and two taking-value options (-s
 * --buffer-size, -n), you can turn on flag a and v, set option s as 100 and n
 * as 50, and pass another two non-option arguments (abc and xyz ) to foo by
 * these commands:\n
 * <tt>
 * \$foo -av -s 100 -n 50 abc xyz \n
 * \$foo -av --buffer-size=100 -n 50 abc xyz \n
 * \$foo --all --verbose --buffer-size=100 -n 50 abc xyz \n
 * </tt>
 *
 * POSIX conventions:\n
 * http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/V1_chap12.html\n
 * http://www.gnu.org/s/hello/manual/libc/Argument-Syntax.html\n
 *
 * GNU extensions:\n
 * http://www.gnu.org/prep/standards/html_node/Command_002dLine-Interfaces.html\n
 */

#ifndef __ANYARG_H__
#define __ANYARG_H__

#include <string>
#include <stdlib.h>
#include <string.h>

using std::string;


/**
 * Use this class to parse command line options and arguments.
 * Both short and long options are supported.
 */
class anyarg
{
public:
	/// Construct a anyarg object.
	anyarg();

	/// Destructor.
	~anyarg();

	/**
	 * Define a boolean flag for a program.
	 * @param letter Single-letter label of a flag.
	 * @param name   Long name of a flag, using hypen to join multiple wrods.
	 * @param help   Description of a flag, will be used to generate usage.
	 * @return \c true if success, \c false if the flag has been defined.
	 */
	bool set_flag(const char *name, char letter, const char *help);

	/// Define a boolean flag only with a short name.
	bool set_flag(char letter, const char *help);

	/// Define a boolean flag only with a long name.
	bool set_flag(const char *name, const char *help);

	/**
	 * Define a taking-value option for a program.
	 * @param name   Long name of an option, using hyphen to join multiple wrods.
	 * @param letter Single-letter label of an option.
	 * @param vtype  The value type of an option, a c-string in uppercase.
	 * @param help   Description of an option, will be used to generate usage.
	 * @return \c true if success, \c false if the flag has been defined.
	 *
	 * Commonly used vtype: "FILE", "DIR", "STRING", "NUM", "SIZE", "INT", "FLOAT".
	 */
	bool set_option(const char *name, char letter, const char *vtype, const char *help);

	/// Define a taking-value option only with a short name.
	bool set_option(char letter, const char *vtype, const char *help);

	/// Define a taking-value option only with a long name.
	bool set_option(const char *name, const char *vtype, const char *help);

	/**
	 * Parse command line arguments.
	 * argv[0] will be parsed as the name of the program, so you should ensure
	 * that argv[0] is the name of the program. Otherwise, the command line
	 * will be wrongly parsed.
	 * @pre Flags and options have to be defined before parsing the command line.
	 * @param argc The count of arguments on command line.
	 * @param argv Command line arguments as a c-string array.
	 * @return \c true if success, \c false if fail.
	 */
	bool parse_command_line(int argc, char **argv);

	/**
	 * Get the value of a flag.
	 * @param letter Single-letter label of a flag.
	 * @return \c true if a flag is specified on command line, \c false if not.
	 */
	bool get_flag(char letter) const;

	/**
	 * Get the value of a flag.
	 * @param name Long name of a flag.
	 * @return \c true if a flag is specified on command line, \c false if not.
	 */
	bool get_flag(const char *name) const;

	/**
	 * Get the value of a taking-value option.
	 * @param letter Single-letter label of an option.
	 * @return The value of an option if it was specified on command line, NULL if not.
	 */
	const char *get_value(char letter) const;

	/**
	 * Get the value of a taking-value option.
	 * @param name Long name of an option.
	 * @return The value of an option if it was specified on command line, NULL if not.
	 */
	const char *get_value(const char *name) const;

	/**
	 * Get the count of non-option arguments.
	 * return The count of non-option arguments.
	 */
	int get_argc() const;

	/**
	 * Get a non-option argument.
	 * @param i The index of a non-option arguments, ranging from 0 to get_argc() - 1.
	 * Non-option arguments are in the order as they appeared on command line.
	 */
	const char *get_arg(int i) const;

	/**
	 * Generate formatted help for program options.
	 * @return Help information as a c-string.
	 */
	const char *autohelp();

	/// Display program options and values, for debugging.
	void display() const;

private:
	char *prog_;  // name of the program
	int argc_;    // count of non-option arguments
	char **argv_; // array of non-option arguments

	// structure of an option
	struct option {
		char type;    // 1 for binary flag, 2 for option with a value
		char letter;  // single-letter label of an option
		char *name;   // long name of an option
		char *value;  // value of an option
		char *vtype;  // value type of an option, such as FILE, NUM, INT, SIZE, in uppercase
		char *help;   // help information about an option
	};

	option *options_; // option array of a program
	int n_option_;    // number of options
	int m_option_;    // the maximum of options that anyarg::options_ can store

	char prefix_letter_;  // prefix of single-letter option
	char prefix_name_[2]; // prefix of long option

	string help_; // formatted help of a program

	anyarg(const anyarg &); // prevent copy

	anyarg & operator = (const anyarg &); // prevent assignment

	bool add_option(int type, const char *name, char letter, const char *vtype, const char *help);

	int get_optind(char letter) const;

	int get_optind(const char *name) const;

	bool set_value(option &option, const char *value);
};

#endif

/**
\example
\code
#include <stdio.h>
#include "anyarg.h"

int main (int argc, char **argv)
{
	// step 1: create a anyarg object.
	anyarg opt;

	// step 2: define usable options
	opt.set_flag("all", 'a', "Include everything");
	opt.set_flag("help", 'h', "Give help information");
	opt.set_flag('v', "Open the verbose mode");
	opt.set_option('c', "COUNT", "Count of books");
	opt.set_option("nind", 'n', "NUM", "Number of students");
	opt.set_option("size", 's', "INT", "Size of a basketball asketball");
	opt.set_option("width", "SIZE", "Width of a picture");
	opt.set_option("long-long-long", "UNKNOWN", "no description");

	// step 3: parse command line arguments
	bool ret;
	ret = opt.parse_command_line(argc, argv);
	if (ret == false) {
		fprintf(stderr, "Fail to parse command line arguments\n");
		return 1;
	}


	// debug, display option and option value
//	opt.display();

	// step 4: get option values and non-option arguments
	const char *s;
	if (opt.get_flag('h')) printf("h is turned on\n");
	if (opt.get_flag('v')) printf("v is turned on\n");
	if (opt.get_flag("all")) printf("all is turned on\n");
	if ((s = opt.get_value("size")) != NULL) printf("size is %d\n", atoi(s));

	int k = opt.get_argc();
	printf("%d non-option arguments\n", k);
	for (int i = 0; i < k; i++)
		printf("arg %d is: %s\n", i, opt.get_arg(i));

	// step 5: get formatted option help if required.
	if (opt.get_flag('h')) printf("%s\n", opt.autohelp());

	return 0;
}
\endcode
*/


