#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "anyarg.h"

using std::string;

#define ANYARG_TYPE_FLAG 1
#define ANYARG_TYPE_OPT  2

// remove redundant spaces in a sentence.
static void remove_space(char *s)
{
	if (s == NULL || s[0] == 0) return;
	int slen = strlen(s);
	int start = 0, end, w;
	int dest = 0;
	int is_first_word = 1;
	while (start < slen) {
		// get a new word
		while (start < slen && isspace(s[start])) start++; // skip leading space
		if (start >= slen) break;
		end = start + 1;
		while (end < slen && isspace(s[end]) == 0) end++; // goto word end
		w = end - start;

		if (is_first_word) {
			is_first_word = 0;
		} else s[dest++] = ' ';

		if (dest != start) memmove(s + dest, s + start, w);
		dest += w;
		start = end;
	}
	s[dest] = '\0';
}

anyarg::anyarg()
{
	prog_ = NULL;
	argc_ = 0; argv_ = NULL;
	options_ = NULL; n_option_ = 0;
	prefix_letter_ = '-';
	prefix_name_[0] = prefix_name_[1] = '-';
	help_ = "";

	m_option_ = 16;
	options_ = (option *) calloc(m_option_, sizeof(option));
	if (options_ == NULL) {
		fprintf(stderr, "%s line %d: fail to allocate memory\n", __FILE__, __LINE__);
		abort();
	}
}

anyarg::~anyarg()
{
	for (int i = 0; i < n_option_; i++) {
		if (options_[i].name)  free(options_[i].name);
		if (options_[i].value) free(options_[i].value);
		if (options_[i].vtype) free(options_[i].vtype);
		if (options_[i].help)  free(options_[i].help);
	}
	free(options_);

	for (int i = 0; i < argc_; i++) free(argv_[i]);
	free(argv_);

	free(prog_);
}

bool anyarg::add_option(int type, const char *name, char letter, const char *vtype, const char *help)
{
	if ((name == NULL || name[0] == 0) && letter == 0) {
		fprintf(stderr, "%s line %d: invalid option\n", __FILE__, __LINE__);
		return 0;
	}

	// check whether the option has been defined.
	int found = 0;
	if (name && name[0]) {
		int j = get_optind(name);
		if (j >= 0) found = 1;
	}
	if (letter) {
		int j = get_optind(name);
		if (j >= 0) found = 1;
	}
	if (found) {
		fprintf(stderr, "%s line %d: option has been defined\n", __FILE__, __LINE__);
		return 0;
	}

	// add option
	if (n_option_ >= m_option_) {
		m_option_ += 16;
		options_ = (option *) realloc(options_, m_option_ * sizeof(option));
	}
	option &t = options_[n_option_];
	t.type = type;
	t.letter = letter;
	if (name && name[0]) t.name = strdup(name);
	else t.name = NULL;
	if (vtype && vtype[0]) t.vtype = strdup(vtype);
	else t.vtype = strdup("UNKNOWN");
	if (help && help[0]) t.help = strdup(help);
	else t.help = strdup("No description avaliable");
	remove_space(t.help);
	t.value = NULL;
	n_option_++;
	return 1;
}

bool anyarg::set_flag(const char *name, char letter, const char *help) {return add_option(ANYARG_TYPE_FLAG, name, letter, "BOOL", help);}
bool anyarg::set_flag(char letter, const char *help) {return add_option(ANYARG_TYPE_FLAG, NULL, letter, "BOOL", help);}
bool anyarg::set_flag(const char *name, const char *help) {return add_option(ANYARG_TYPE_FLAG, name, 0, "BOOL", help);}

bool anyarg::set_option(const char *name, char letter, const char *vtype, const char *help) {return add_option(ANYARG_TYPE_OPT, name, letter, vtype, help);}
bool anyarg::set_option(char letter, const char *vtype, const char *help) {return add_option(ANYARG_TYPE_OPT, NULL, letter, vtype, help);}
bool anyarg::set_option(const char *name, const char *vtype, const char *help) {return add_option(ANYARG_TYPE_OPT, name, 0, vtype, help);}

// Get the index of a simple option
// return -1 if not found
int anyarg::get_optind(char letter) const
{
	if (letter == 0) return -1;
	for (int i = 0; i < n_option_; i++)
		if (letter == options_[i].letter) return i;
	return -1;
}

// Get the index of a long option
int anyarg::get_optind(const char *name) const
{
	if (name == NULL || name[0] == 0) {
		return -1;
	}
	for (int i = 0; i < n_option_; i++) {
		if (options_[i].name == NULL) continue;
		if (strcmp(name, options_[i].name) == 0) return i;
	}
	return -1;
}

// Turn on flags specified on command line
// Get option values specified on command line
// Collect non-option arguments
bool anyarg::parse_command_line(int argc, char **argv)
{
	prog_ = strdup(argv[0]);

	argc_ = 0;
	argv_ = (char **) calloc(argc, sizeof(char *));

	int stop_parsing = 0; // stop option parsing when meeting argument --
	int i, j;
	for (i = 1; i < argc; i++) {
		if (stop_parsing) {argv_[argc_++] = strdup(argv[i]); continue;}
		if (argv[i][0] == prefix_name_[0] && argv[i][1] == prefix_name_[1]) { // long option
			if (argv[i][2] == 0) {stop_parsing = 1; continue;}
			char *o = strdup(argv[i] + 2);
			int olen = strlen(o);
			int t = 0;
			while (o[t] && o[t] != '=') t++; // seek for equal sign
			if (t < olen) { // is a option=value pair
				o[t] = '\0';
				j = get_optind(o);
				if (j < 0) { // the option is undefined
					fprintf(stderr, "%s line %d: long option is undefined\n", __FILE__, __LINE__);
					return 0;
				}
				if (options_[j].type == ANYARG_TYPE_FLAG) {
					fprintf(stderr, "%s line %d: cannot give a value to a flag\n", __FILE__, __LINE__);
					return 0;
				}
				set_value(options_[j], o + t + 1);
			} else { // is a regular option
				j = get_optind(o);
				if (j < 0) {
					fprintf(stderr, "%s line %d: long option is undefined\n", __FILE__, __LINE__);
					return 0;
				}
				if (options_[j].type == ANYARG_TYPE_FLAG) {
					set_value(options_[j], "ON");
				} else {
					if (i >= argc - 1) {
						fprintf(stderr, "%s line %d: option value does not exist\n", __FILE__, __LINE__);
						return 0;
					}
					set_value(options_[j], argv[++i]);
				}
			}
			free(o);
		} else if (argv[i][0] == prefix_letter_) { // short option
			if (argv[i][1] == 0) {argv_[argc_++] = strdup("-"); continue;}
			const char *o = argv[i] + 1;
			int olen = strlen(o);
			if (olen == 1) { // is a regular short option
				j = get_optind(o[0]);
				if (j < 0) {
					fprintf(stderr, "%s line %d: short option is undefined\n", __FILE__, __LINE__);
					return 0;
				}
				if (options_[j].type == ANYARG_TYPE_FLAG) {
					set_value(options_[j], "ON");
				} else {
					if (i >= argc - 1) {
						fprintf(stderr, "%s line %d: option value does not exist\n", __FILE__, __LINE__);
						return 0;
					}
					set_value(options_[j], argv[++i]);
				}
			} else { // bundled flags ('-lst') or option-value pair ('-ofoo');
				/// no bundled flags, because this design will conflict with wrong long options.
				fprintf( stderr, "%s long option is undefined\n", argv[i] );
				abort();
				j = get_optind(o[0]);
				if (j < 0) {
					fprintf(stderr, "%s line %d: short option is undefined\n", __FILE__, __LINE__);
					return 0;
				}
				if (options_[j].type == ANYARG_TYPE_OPT) { // is option-value pair
					set_value(options_[j], o+1);
				} else { // is bundled flags
					for (int t = 0; t < olen; t++) {
						j = get_optind(o[t]);
						if (j < 0) {
							fprintf(stderr, "%s line %d: short option is undefined\n", __FILE__, __LINE__);
							return 0;
						}
						if (options_[j].type != ANYARG_TYPE_FLAG) {
							fprintf(stderr, "%s line %d: only short flags can be grouped\n", __FILE__, __LINE__);
							return 0;
						}
						set_value(options_[j], "ON");
					}
				}
			}
		} else { // non-option arguments
			argv_[argc_++] = strdup(argv[i]);
		}
	}
	return 1;
}

bool anyarg::set_value(option &opt, const char *val)
{
	if (val == NULL) {
		fprintf(stderr, "%s line %d: option value cannot be NULL\n", __FILE__, __LINE__);
		return 0;
	}
	if (opt.value != NULL) {
		fprintf(stderr, "%s line %d: overwritten of option value is disallowed\n", __FILE__, __LINE__);
		return 0;
	}
	opt.value = strdup(val);
	return 1;
}

bool anyarg::get_flag(char letter) const
{
	int j = get_optind(letter);
	if (j < 0) {
		fprintf(stderr, "%s line %d: short option is undefined\n", __FILE__, __LINE__);
		return 0;
	}
	if (options_[j].type != ANYARG_TYPE_FLAG) {
		fprintf(stderr, "%s line %d: not a flag\n", __FILE__, __LINE__);
		return 0;
	}
	if (options_[j].value && strcmp(options_[j].value, "ON") == 0) return 1;
	return 0;
}

const char *anyarg::get_value(char letter) const
{
	int j = get_optind(letter);
	if (j < 0) {
		fprintf(stderr, "%s line %d: short option is undefined\n", __FILE__, __LINE__);
		return 0;
	}
	if (options_[j].type != ANYARG_TYPE_OPT) {
		fprintf(stderr, "%s line %d: not a option\n", __FILE__, __LINE__);
		return 0;
	}
	return options_[j].value;
}

bool anyarg::get_flag(const char *name) const
{
	int j = get_optind(name);
	if (j < 0) {
		fprintf(stderr, "%s line %d: long option is undefined\n", __FILE__, __LINE__);
		return 0;
	}
	if (options_[j].type != ANYARG_TYPE_FLAG) {
		fprintf(stderr, "%s line %d: not a flag\n", __FILE__, __LINE__);
		return 0;
	}
	if (options_[j].value && strcmp(options_[j].value, "ON") == 0) return 1;
	return 0;
}

const char *anyarg::get_value(const char *name) const
{
	int j = get_optind(name);
	if (j < 0) {
		fprintf(stderr, "%s line %d: long option is undefined\n", __FILE__, __LINE__);
		return NULL;
	}
	if (options_[j].type != ANYARG_TYPE_OPT) {
		fprintf(stderr, "%s line %d: not a option\n", __FILE__, __LINE__);
		return NULL;
	}
	return options_[j].value;
}

int anyarg::get_argc() const {return argc_;}

const char *anyarg::get_arg(int i) const
{
	if (i < 0 || i >= argc_) {
		fprintf(stderr, "%s line %d: index out of range\n", __FILE__, __LINE__);
		return NULL;
	}
	return argv_[i];
}

// fix a bug here
const char *anyarg::autohelp()
{
	int nindent = 24;

	help_.clear();
	help_.append("Options\n");

	for (int i = 0; i < n_option_; i++) {
		option &opt = options_[i];
		int k = 0;
		char s[256]; // keep formatted option string
		for (int j = 0; j < 256; j++) s[j] = ' ';
		if (opt.letter) {
			s[2] = '-'; s[3] = opt.letter;
		}
		if (opt.name && opt.name[0]) {
			if (opt.letter) s[4] = ',';
			s[6] = s[7] ='-';
			int l = strlen(opt.name);
			if (l > 64) l = 64; // for safe, impossible to happen
			k = 8;
			strncpy(s + k, opt.name, l); k += l;
			if (opt.type == ANYARG_TYPE_OPT) {
				s[k++] = '=';
				l = strlen(opt.vtype);
				if (l > 64) l = 64; // for safe, impossible to happen
				strncpy(s+k, opt.vtype, l); k += l;
			}
		}
		if (k < nindent - 2) {
			s[nindent] = '\0';
			help_.append(s);
		} else {
			s[k] = '\0';
			help_.append(s);
			help_.append("\n");
			help_.append(nindent, ' ');
		}
		help_.append(opt.help);
		help_.append("\n");
	}
	return help_.c_str();
}

void anyarg::display() const
{
	int i;
	printf("options\n");
	printf("letter | name | value | vtype | help\n");
	for (i = 0; i < n_option_; i++) {
		option &o = options_[i];
		printf("%c\t%s\t%s\t%s\t%s\n", o.letter, o.name, o.value, o.vtype, o.help);
	}
	printf("\n");

	printf("\n%d non-option arguments\n", argc_);
	for (i = 0; i < argc_; i++) printf("  %s\n", argv_[i]);
}

