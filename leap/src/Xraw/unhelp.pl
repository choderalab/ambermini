: # use perl                                  -*- mode: Perl; -*-
        eval 'exec perl -S $0 "$@"'
                if $runnning_under_some_shell;

foreach $file (@ARGV) {
    open(FILE, $file) || die ("Cannot open file $file");
    $save = "./new/" . $file;

    open(SAVE, "> $save") || die ("Cannot open file $save");

    @contents = <FILE>;
    
    print "making $save..\n";

    foreach (@contents) {
	if (! /^\#include <dmalloc.h>/) {
	    print SAVE $_;
	}
    }
	
    close(FILE);
    close(SAVE);
}
