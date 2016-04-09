// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once
#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <assert.h>
#include <math.h>
#include <direct.h>
#include <time.h>

#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>

//生成一个[0,1)的实数
#define RANDI (rand()/(RAND_MAX+1.0))
//生成一个[0,MAX)的整数
#define RANDIUINT(MAX) ((int)floor(RANDI * (MAX)))

#ifdef _DUBUG
#define LOG(...) printf(__VA_ARGS__)
#else
#define LOG(...)
#endif

#define CFree(_Memory) free(_Memory); _Memory = NULL;
#define CDelete(_Memory) delete _Memory; _Memory = NULL;
#define CDeleteArray(_Memory) delete[] _Memory; _Memory =NULL;
//根据极坐标转换成直角坐标X,Y,Z
//http://baike.baidu.com/link?url=zbluOmziAKGRCfxSU6DwjJ3dSp2wyTzLLGBJEZpN8yUE4vCTk7fPnAa0C2KMcL81ImoROM24R2N4OVYq5n2BSa
#define XFROMPOLE(r,theta,phi) (r*sin(theta)*cos(phi))
#define YFROMPOLE(r,theta,phi) (r*sin(theta)*sin(theta))
#define ZFROMPOLE(r,theta,phi) (r*cos(theta))

#ifndef M_PI
	#define M_PI       3.14159265358979323846
#endif

using namespace std;


// TODO: 在此处引用程序需要的其他头文件
