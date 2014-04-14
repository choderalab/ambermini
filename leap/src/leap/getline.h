#ifndef GETLINE_H
#define GETLINE_H

/* unix systems can #define POSIX to use termios, otherwise 
 * the bsd or sysv interface will be used 
 */

extern char	*tl_getline(char *prompt);		/* read a line of input */
void            gl_histadd(char *buf);		/* adds entries to hist */

#endif /* GETLINE_H */
