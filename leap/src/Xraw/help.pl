: # use perl                                  -*- mode: Perl; -*-
        eval 'exec perl -S $0 "$@"'
                if $runnning_under_some_shell;

foreach $file (@ARGV) {
    open(FILE, $file) || die ("Cannot open file $file");
    $save = "./new/" . $file;

    open(SAVE, "> $save") || die ("Cannot open file $save");

    @contents = <FILE>;
    
    print "making $save..\n";

    $line = 0;
    for($i = $#contents; $i > $[; --$i) {
	if ($contents[$i] =~ /^\#include/) {
	    $line = $i;
	    last;
	}
    }
	
    for ($j = $[; $j <= $line; ++$j) {
	print SAVE $contents[$j];
    }
    
    print SAVE "\n\#include <dmalloc.h>\n";
    
    for ($j = $line+1; $j <= $#contents; ++$j) {
	print SAVE $contents[$j];
    }

    close(FILE);
    close(SAVE);
}
