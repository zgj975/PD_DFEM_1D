
// PD Manager.h : PROJECT_NAME Ӧ�ó������ͷ�ļ�
//

#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"		// ������


// CPDManagerApp: 
// �йش����ʵ�֣������ PD Manager.cpp
//

class CPDManagerApp : public CWinApp
{
public:
	CPDManagerApp();

// ��д
public:
	virtual BOOL InitInstance();
	virtual BOOL ExitInstance();
// ʵ��

	DECLARE_MESSAGE_MAP()
#ifdef _DEBUG
protected:
	CMemoryState m_msOld, m_msNew, m_msDiff;
#endif // _DEBUG
};

extern CPDManagerApp theApp;