#ifndef MyProcessManager_h
#define MyProcessManager_h 1

#include <stdio.h>

class MyProcessManager{
	public:
		MyProcessManager();
		virtual ~MyProcessManager();

		static MyProcessManager* GetMyProcessManager();

		int OpenFile();
		void CloseFile();

		double GetMemorySize();

	private:
		static MyProcessManager* fMyProcessManager;
		int MemoryConsumption;
		FILE *fin_card;
		bool fileOpened;
};

#endif
