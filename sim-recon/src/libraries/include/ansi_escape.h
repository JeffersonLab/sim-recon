
// These are ANSI escape codes used to modify how characters
// are drawn on a vt100-like terminal (e.g. xterm).
//
// Example 1:
//    Write a string in bold font
//
//     cout<<ansi_bold<<"I am bold!"<<ansi_normal<<endl;
//
//
// Example 2:
//    Write a message in green and then go up 2 lines and
//    write another in red
//
//     cout<<ansi_green<<"Spring is wonderful"<<ansi_normal<<endl;
//     cout<<ansi_up(3)<<ansi_red<<"We're all going to die!"<<ansi_normal<<ansi_down(3);
//
//
//  Note that the character drawing attributes (bold, color, etc...)
// will stay with the terminal even after your program exits. As such,
// it is advisable to immediately send a ansi_normal sequence to the 
// screen as soon as you've finished writing  the formatted message.
//

#ifndef ansi_escape
#define ansi_escape			((char)0x1b)
#define ansi_bold 			ansi_escape<<"[1m"
#define ansi_italic 			ansi_escape<<"[3m"
#define ansi_underline 		ansi_escape<<"[4m"
#define ansi_blink 			ansi_escape<<"[5m"
#define ansi_rapid_blink	ansi_escape<<"[6m"
#define ansi_reverse			ansi_escape<<"[7m"
#define ansi_black			ansi_escape<<"[30m"
#define ansi_red				ansi_escape<<"[31m"
#define ansi_green			ansi_escape<<"[32m"
#define ansi_blue				ansi_escape<<"[34m"
#define ansi_normal			ansi_escape<<"[0m"
#define ansi_up(A)			ansi_escape<<"["<<(A)<<"A"
#define ansi_down(A)			ansi_escape<<"["<<(A)<<"B"
#define ansi_forward(A)		ansi_escape<<"["<<(A)<<"C"
#define ansi_back(A)			ansi_escape<<"["<<(A)<<"D"
#endif // ansi_escape
