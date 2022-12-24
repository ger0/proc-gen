#pragma once

// only temporary
#define DEBUG

// logs
#ifdef DEBUG
constexpr const char* NRM = "\x1B[0m";
constexpr const char* RED = "\x1B[31m";
constexpr const char* GRN = "\x1B[32m";
constexpr const char* YEL = "\x1B[33m";
constexpr const char* BLU = "\x1B[34m";
constexpr const char* MAG = "\x1B[35m";
constexpr const char* CYN = "\x1B[36m";
constexpr const char* WHT = "\x1B[37m";

#define LOG_INFO(STR, ...) 	printf("%s" STR "%s\n", YEL, ##__VA_ARGS__, WHT)
#define LOG_PRINT(STR) 		printf("%s%s%s\n", MAG, STR, WHT)
#define LOG_WARN(STR, ...) 	fprintf(stderr, "%s" STR "%s\n", RED, ##__VA_ARGS__, WHT)

#else 
#define LOG_INFO(...) ;
#define LOG_WARN(...) ;
#define LOG_PRINT(...) ;
#endif
