//#**************************************************************
//#
//# filename:             TcpIpRoutines.cpp
//#
//# author:               Schörgenhumer Markus
//#					
//# generated:						Feb and July 2013
//#
//# description:          - base class TCPIPSocket for TCP/IP server/client
//#													socket set-up, communication, and error handling
//#												- derived class TCPIPHotInt for use in HOTINT 
//#
//# remarks:              - additional linking dependencies: ws2_32.lib
//# 											  (in HOTINT in the project "MBSElementsAndModels")
//#
//# Copyright (c) 2003-2013 Johannes Gerstmayr, Linz Center of Mechatronics GmbH, Austrian
//# Center of Competence in Mechatronics GmbH, Institute of Technical Mechanics at the 
//# Johannes Kepler Universitaet Linz, Austria. All rights reserved.
//#
//# This file is part of HotInt.
//# HotInt is free software: you can redistribute it and/or modify it under the terms of 
//# the HOTINT license. See folder 'licenses' for more details.
//#
//# bug reports are welcome!!!
//# WWW:		www.hotint.org
//# email:	bug_reports@hotint.org or support@hotint.org
//#**************************************************************


#include <fstream>
#include <iostream>
#include <sstream>
#include "TcpIpRoutines.h"
#include "winsock2.h"        //note: also includes windows.h

using namespace std;


//========================================================================================
//================================== TCPIPSocket =========================================
//========================================================================================

TCPIPSocket::TCPIPSocket(SocketType typei, char* ipdata, int errormsgflag, int infomsgflag)
{
	type=typei;
	port=0;
	//ipstring=NULL;
	status=0;
	autoreconnect=0;
	looptime=1000;
	//recvbuffer=0;
	SetErrorMessageMode(errormsgflag);
	SetInfoMessageMode(infomsgflag);
	GetTCPIPConfig(ipdata);
	accepttimeout = -1;
}

TCPIPSocket::TCPIPSocket(SocketType typei, char* ipIn, unsigned short portIn, int errormsgflag, int infomsgflag)
{
	Set(typei,ipIn,portIn,errormsgflag,infomsgflag);
}

void TCPIPSocket::Set(SocketType typei, char* ipIn, unsigned short portIn, int errormsgflag, int infomsgflag)
{
	type=typei;
	port=portIn;
	//ipstring=ipIn;
	strcpy(ipstring,ipIn);
	status=0;
	autoreconnect=0;
	looptime=1000;
	SetErrorMessageMode(errormsgflag);
	SetInfoMessageMode(infomsgflag);
	accepttimeout = -1;
}


int TCPIPSocket::SetUpConnection()
{	
	if(status)
		CloseConnection();

	if(type==ST_client)
	{
		StartWinsock();
		CreateSocket();
		//ConnectOnce();
		ConnectLoop();
		//FlushRecvBuffer();
		Initialize();
		return 1;
	}
	else if(type==ST_server)
	{
		StartWinsock();
		CreateSocket();
		BindSocket();
		ListenMode();
		AcceptConnection();
		Initialize();
		return 1;
	}
	else
	{
		Error(0,"Invalid socket type!");
		return -1;
	}
}

void TCPIPSocket::CloseConnection()
{
	if(status)
	{
		closesocket(c);
		if(type==ST_server)
			closesocket(s);
		WSACleanup();
		status = 0;
	}
}

TCPIPSocket::~TCPIPSocket()
{
	CloseConnection();
	//if(ipstring)
	//	delete [] ipstring;
}


void TCPIPSocket::GetTCPIPConfig(char* file)
{
	ifstream ipconfig;
	string a1,a2,a3,a4,temp1;

	ipconfig.open(file,ios::binary|ios::in);
	if(!ipconfig.is_open())
	{
		Error(0,string("Access to TCP/IP configuration file '").append(string(file)).append(string("' failed!")));
		//cout << "Error accessing  'IP.txt'" << endl;
		//cout << "Enter configuration manually: " << endl << endl;
		//cout << "Server-IP a1.a2.a3.a4 (IPv4): " << endl << endl;
		//cout << "a1: ";
		//cin >> a1;
		//cout << "a2: ";
		//cin >> a2;
		//cout << "a3: ";
		//cin >> a3;
		//cout << "a4: ";
		//cin >> a4;
		//cout << "Port: ";
		//cin >> port;
		//cout << endl;
	}
	else
	{
		ipconfig >> a1;
		ipconfig >> a2;
		ipconfig >> a3;
		ipconfig >> a4;
		ipconfig >> port;
		ipconfig.clear();
		ipconfig.close();
	}

	temp1.append(a1).append(".").append(a2).append(".").append(a3).append(".").append(a4);
	//ipstring = new char[temp1.length()+1];   //generation of a c-strings for the IP
	strcpy(ipstring,temp1.c_str());

	//stringstream out;
	//out << "server-IP / port: " << ipstring << " / " << port << endl << endl;
	//InfoMessage(out.str());
}

long TCPIPSocket::StartWinsock(){
  WSADATA wsa;
  long rc = WSAStartup(MAKEWORD(2,0),&wsa);  //MAKEWORD ... macro which transforms the version number 2.0 into WORD (unsigned short)
	if(rc)
		Error(1,"Failed to start winsock!",rc);
	else
		InfoMessage("Winsock started");
	return rc;
}

void TCPIPSocket::CreateSocket()
{
	//SOCKET & loc = s;    //somehow, this reference IS NOT WORKING, no clue why; workaround below the commented lines
	//if(type==ST_client)
	//	loc=c;
	//loc=socket(AF_INET,SOCK_STREAM,0);
	//if(loc==INVALID_SOCKET)
	//	Error(1,"Failed to create socket!\n",WSAGetLastError());
	//else
	//	InfoMessage("Socket created\n");

	if(type==ST_client)
	{
		c=socket(AF_INET,SOCK_STREAM,0);
		if(c==INVALID_SOCKET)
			Error(1,"Failed to create socket!",WSAGetLastError());
		else
			InfoMessage("Socket created");
	}
	else
	{
		s=socket(AF_INET,SOCK_STREAM,0);
		if(s==INVALID_SOCKET)
			Error(1,"Failed to create socket!",WSAGetLastError());
		else
			InfoMessage("Socket created");
	}
}

long TCPIPSocket::BindSocket()
{
	SOCKADDR_IN addr;
	long rc;

	memset(&addr,0,sizeof(SOCKADDR_IN));   //cf. p.19, Beej's guide
	addr.sin_family=AF_INET; //IPv4
	addr.sin_port=htons(port);
	addr.sin_addr.s_addr=inet_addr(ipstring);
	rc=bind(s,(SOCKADDR*)&addr,sizeof(SOCKADDR_IN));
	if(rc==SOCKET_ERROR)
		Error(1,"Failed to bind server socket to given IP and port!",WSAGetLastError());
	else
	{
		stringstream out;
		out << "Server socket bound to IP " << ipstring << ", port " << port;
		InfoMessage(out.str());
	}
	return rc;
}

long TCPIPSocket::ListenMode()
{
	long rc;
	rc=listen(s,10);
	if(rc==SOCKET_ERROR)
		Error(1,"Failed to activate listen mode!",WSAGetLastError());
	else
		InfoMessage("Server socket is in listen mode...");
	return rc;
}

int TCPIPSocket::AcceptConnection()
{
	if(accepttimeout>0)
	{
		int i=0;
		while(1)
		{
			//workaround to implement timeout for accept

			//from MSDN on the "select"-function:
			//The parameter readfds identifies the sockets that are to be checked for readability. If the socket is currently in the listen state, it will be marked as readable
			//if an incoming connection request has been received such that an accept is guaranteed to complete without blocking.

			timeval to;
			to.tv_sec = (unsigned int) ((double)accepttimeout)/1000;
			to.tv_usec = accepttimeout%1000;
			fd_set socketset;
			socketset.fd_count = 1;
			socketset.fd_array[0] = s;
			//FD_ZERO(&socketset);
			//FD_SET(c,&socketset);
			int temp = select(1,&socketset,NULL,NULL,&to);
			if(temp==0)
			{
				++i;
				if(i<3)
					Message("Failed to accept connection from client! (timeout) Retry...");
				else
					Error(1,"Failed to accept connection from client! (timeout)");
			}
			else if(temp==SOCKET_ERROR)
				Error(1,"Failed to accept connection from client!",WSAGetLastError());
			else
				break;
		}
	}

	c=accept(s,NULL,NULL);
	if(c==INVALID_SOCKET)
		Error(1,"Failed to accept connection from client!",WSAGetLastError());
	else
	{
		InfoMessage("New connection to client accepted");
		status=1;
	}
	return 1;
}

long TCPIPSocket::ConnectLoop()
{
	long rc;
	SOCKADDR_IN addr;	//addr contains all neccessary information about the server socket; the port of the client socket is chosen by Windows
	stringstream out;

	memset(&addr,0,sizeof(SOCKADDR_IN)); //reset everything to 0
	addr.sin_family=AF_INET;
	addr.sin_port=htons(port); //port
	addr.sin_addr.s_addr=inet_addr(ipstring); //IP v4 adress of server
	
	int waitingcounter = -1;
	while(1){
		rc=connect(c,(SOCKADDR*)&addr,sizeof(SOCKADDR));
		if(rc==SOCKET_ERROR)
		{
			switch(waitingcounter){
					case -1: out.str(""); out << "Waiting to connect to " << ipstring << ", port " << port << endl; break;
					case 0:  out.str(""); out << "\r." << "   "; break;
					case 1:  out.str(""); out << "\r.." << "  "; break;
					case 2:  out.str(""); out << "\r..." << " "; waitingcounter = -1; break;
			}
			InfoMessage(out.str());
			++waitingcounter;
			Sleep(looptime); //wait for "time" milliseconds 
		}
		else
		{
			InfoMessage("Connected!");
			status=1;
			break;
		}
	}
	return rc;
}

long TCPIPSocket::ConnectOnce()
{
	long rc;
	SOCKADDR_IN addr;	//addr contains all neccessary information about the server socket; the port of the client socket is chosen by Windows
	
	memset(&addr,0,sizeof(SOCKADDR_IN)); //reset everything to 0
	addr.sin_family=AF_INET;
	addr.sin_port=htons(port); // port
	addr.sin_addr.s_addr=inet_addr(ipstring); //IP v4 adress of server
	rc=connect(c,(SOCKADDR*)&addr,sizeof(SOCKADDR));

	if(rc==SOCKET_ERROR)
		Error(1,"Connect failed!",WSAGetLastError());
	else
	{
		InfoMessage("Connected!");
		status=1;
	}
	return rc;
}

//void TCPIPSocket::FlushRecvBuffer()
//{
//	int bufsize;
//	int optlen=sizeof(int);
//	getsockopt(c,sol_socket,so_rcvbuf,reinterpret_cast<char*> (&bufsize),&optlen);
//	//cout << "recv buffer in bytes: " << temp << "\n";
//	char* temp = new char[bufsize];
//
//	//unsigned long recvto;
//	//optlen=sizeof(unsigned long);
//	//getsockopt(c,sol_socket,so_rcvtimeo,reinterpret_cast<char*> (&recvto),&optlen);
//	//cout << "timeout " << recvto << "\n";
//
//	setrecvtimeout(10000);
//
//	int remain=bufsize;
//	int tmp=0;
//	int count=0;
//	while(remain!=0 && tmp>=0)
//	{
//		int tmp = recv(c,temp+(bufsize-remain),remain,0);
//		remain-=tmp;
//		count+=tmp;
//	}
//
//	cout << "from buffer: " << count << "\n";
//	cout << "last error" << wsagetlasterror() << "\n";
//
//	//setrecvtimeout(recvto);
//	
//	delete [] temp;
//}

void TCPIPSocket::Initialize()
{
}

int TCPIPSocket::SendByte(byte data)
{
	int stat = SendData(reinterpret_cast<char*> (&data),1);
	return stat;
}

int TCPIPSocket::RecvByte(byte & data)
{
	int stat = RecvData(reinterpret_cast<char*> (&data),1);	
	return stat;
}

int TCPIPSocket::SendInteger(int data)
{
	int temp=htonl(data);
	int stat = SendData(reinterpret_cast<char*> (&temp),4);
	return stat;
}

int TCPIPSocket::RecvInteger(int & datain)
{
	datain=0;
	int stat = RecvData(reinterpret_cast<char*> (&datain),4);	
	datain=ntohl(datain);
	return stat;
}

int TCPIPSocket::SendData(char* data, int size)
{
	if(!status)
	{
		Error(1,"Failed to send data. Socket is not properly initialized!");
		return -1;
	}

	int remain=size;
	int tmp=0;

	while(remain!=0)
	{
		int tmp = send(c,data+(size-remain),remain,0);
		if(tmp<=0)
		{
			Error(1,"An error occurred while sending data!",WSAGetLastError());
			return tmp;
		}
		remain-=tmp;
	}
	return size;
}

int TCPIPSocket::RecvData(char* data, int size)
{
	if(!status)
	{
		Error(1,"Failed to receive data. Socket is not properly initialized!");
		return -1;
	}

	int remain=size;
	int tmp=0;

	while(remain!=0)
	{
		int tmp = recv(c,data+(size-remain),remain,0);
		//if(tmp==0)
		//{
		//	Error(1,"Connection was closed by the server!");
		//	return tmp;
		//}
		if(tmp<=0)
		{
			Error(1,"An error occurred while receiving data!",WSAGetLastError());
			return tmp;
		}
		remain-=tmp;
	}
	return size;
}

void TCPIPSocket::SetRecvTimeOut(unsigned long time)
{
	setsockopt(c,SOL_SOCKET,SO_RCVTIMEO,reinterpret_cast<const char*> (&time),sizeof(unsigned long));
}

void TCPIPSocket::SetSendTimeOut(unsigned long time)
{
	setsockopt(c,SOL_SOCKET,SO_SNDTIMEO,reinterpret_cast<const char*> (&time),sizeof(unsigned long));
}

void TCPIPSocket::SetAcceptTimeOut(unsigned long time)
{
	accepttimeout = time;
}

//void TCPIPSocket::SetRecvBufferSize(int size)
//{
//	recvbuffer=size;
//	setsockopt(c,SOL_SOCKET,SO_RCVBUF,reinterpret_cast<const char*> (&recvbuffer),sizeof(int));
//}

void TCPIPSocket::SetLoopTime(int time)
{
	looptime=time;
}

int TCPIPSocket::Error(int type, string err, int errornumber)
{
	stringstream ts;
	if(type==0)
	{
		if(errormsg)
		{
			string temp = "Fatal error: ";
			if(errornumber)
			{
				ts<<errornumber;
				temp.append(ts.str());
			}
			temp.append("\n");
			std::cout << temp;
			std::cout << err << "\n";
			system("PAUSE");
		}
		exit(-1);
	}
	else if(type==1)
	{
		if(errormsg)
		{
			string temp="Error: ";
			if(errornumber)
			{
				ts<<errornumber;
				temp.append(ts.str());
			}
			temp.append("\n");
			std::cout << temp;
			std::cout << err << "\n";
			system("PAUSE");
		}
		CloseConnection();
		if(autoreconnect)
		{
			InfoMessage("Re-initializing TCP/IP connection...");
			SetUpConnection();
		}
	}
	return errornumber;
}

void TCPIPSocket::Message(string message)
{
	cout << message << endl;
}

void TCPIPSocket::InfoMessage(string info)
{
	if(infomsg)
		cout << info << endl;
}

void TCPIPSocket::SetErrorMessageMode(int flag)
{
	if(flag)
		errormsg=1;
	else
		errormsg=0;
}

void TCPIPSocket::SetInfoMessageMode(int flag)
{
	if(flag)
		infomsg=1;
	else
		infomsg=0;
}

void TCPIPSocket::SetAutoReconnect(int flag)
{
	if(flag)
		autoreconnect=1;
	else
		autoreconnect=0;
}


//========================================================================================
//=================================== TCPIPHotInt ========================================
//========================================================================================

void TCPIPHotInt::Set(SocketType typei, char* ipIn, unsigned short portIn, int errormsgflag, int infomsgflag)
{
	TCPIPSocket::Set(typei,ipIn,portIn,errormsgflag,infomsgflag);
}

void TCPIPHotInt::Set(MBS* mbsI, SocketType typei, char* ipIn, unsigned short portIn, int errormsgflag, int infomsgflag)
{
	mbs=mbsI;
	TCPIPSocket::Set(typei,ipIn,portIn,errormsgflag,infomsgflag);
}

int TCPIPHotInt::Error(int type, std::string err, int errornumber)
{
	stringstream ts;
	if(errormsg)
	{
		string temp = "TCP/IP error: ";
		temp.append(err);
		if(errornumber)
		{
		  temp.append(" (");
			ts<<errornumber;
			temp.append(ts.str());
			temp.append(")\n");
		}

		mbs->UO().InstantMessageText(temp.c_str());
		assert(0);
		exit(0);
	}
	return errornumber;
}

void TCPIPHotInt::Message(std::string message)
{
	mbs->UO() << message.c_str() << "\n";
}

void TCPIPHotInt::InfoMessage(std::string info)
{
	mbs->UO() << info.c_str() << "\n";
}
