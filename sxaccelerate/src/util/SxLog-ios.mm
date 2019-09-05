#import <Foundation/Foundation.h>

void sxLogMsg (const char *tag, const char *msg)
{
   NSLog(@"[%s] %s", tag, msg);
}
