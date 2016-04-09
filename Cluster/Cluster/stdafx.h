// stdafx.h : ��׼ϵͳ�����ļ��İ����ļ���
// ���Ǿ���ʹ�õ��������ĵ�
// �ض�����Ŀ�İ����ļ�
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

//����һ��[0,1)��ʵ��
#define RANDI (rand()/(RAND_MAX+1.0))
//����һ��[0,MAX)������
#define RANDIUINT(MAX) ((int)floor(RANDI * (MAX)))

#ifdef _DUBUG
#define LOG(...) printf(__VA_ARGS__)
#else
#define LOG(...)
#endif

#define CFree(_Memory) free(_Memory); _Memory = NULL;
#define CDelete(_Memory) delete _Memory; _Memory = NULL;
#define CDeleteArray(_Memory) delete[] _Memory; _Memory =NULL;
//���ݼ�����ת����ֱ������X,Y,Z
//http://baike.baidu.com/link?url=zbluOmziAKGRCfxSU6DwjJ3dSp2wyTzLLGBJEZpN8yUE4vCTk7fPnAa0C2KMcL81ImoROM24R2N4OVYq5n2BSa
#define XFROMPOLE(r,theta,phi) (r*sin(theta)*cos(phi))
#define YFROMPOLE(r,theta,phi) (r*sin(theta)*sin(theta))
#define ZFROMPOLE(r,theta,phi) (r*cos(theta))

#ifndef M_PI
	#define M_PI       3.14159265358979323846
#endif

using namespace std;


// TODO: �ڴ˴����ó�����Ҫ������ͷ�ļ�
