# notbracket.pl
# print the lines without the < html > tags 
while (<>) {
  if (/</i | />/i | /\"/i){ 
#    print "hello\n";
    } else {
    if (/\:/){ print ; }
   }
 }   
