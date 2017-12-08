#include "MyProcessManager.h"

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "TString.h"

MyProcessManager* MyProcessManager::fMyProcessManager = 0;

MyProcessManager::MyProcessManager()
{
	if (fMyProcessManager){
	    fprintf(stderr,"ERROR: MyProcessManager already instantiated!\n");
	}
	fileOpened = false;
	fin_card = 0;
	MemoryConsumption = 0;
}

MyProcessManager::~MyProcessManager()
{
	CloseFile();
}

MyProcessManager* MyProcessManager::GetMyProcessManager(){
	if ( !fMyProcessManager ){
		fMyProcessManager = new MyProcessManager;
	}
	return fMyProcessManager;
}

int MyProcessManager::OpenFile(){
	CloseFile();
	pid_t id = getpid();
	fin_card = fopen(Form("/proc/%d/status",id),"r");
	if ( !fin_card ){
		fileOpened = false;
		return 1;
	}
	else{
		fileOpened = true;
		return 0;
	}
}

void MyProcessManager::CloseFile(){
	if (fileOpened&&fin_card){
		fclose(fin_card);
	}
	fileOpened = false;
	fin_card = 0;
}

double MyProcessManager::GetMemorySize(){
	if (!fileOpened){
		OpenFile();
	}
    double size = 0;
	if (fin_card&&fileOpened){
        char * temp;
        TString tempS;
        char buff[1024];
        rewind(fin_card);
		while(fgets(buff,1024,fin_card)!=NULL){
		    temp = strtok(buff," \t\n\r");
		    if (temp) tempS = temp;
			if (tempS=="VmSize:"){
                temp = strtok(NULL," \t\n\r");
                if (temp)
                    size = atoi(temp);
                else
                    size = 0;
                tempS = temp;
                temp = strtok(NULL," \t\n\r");
                if (temp)
                    tempS = temp;
                else
                    tempS = "";
				if (tempS == "kB" || tempS == "KB" || tempS == "kb"){
					size *= 1;
				}
				else if (tempS == "mB" || tempS == "MB" || tempS == "mb"){
					size *= 1024;
				}
				else if (tempS == "gB" || tempS == "GB" || tempS == "gb"){
					size *= 1024*1024;
				}
				else if (tempS == "Byte" || tempS == "B" || tempS == "b"){
					size /= 1024;
				}
				break;
			}
		}
	}
	else{
		size = -1;
	}
	CloseFile();
	return size;
}
