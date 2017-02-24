//#**************************************************************
//#
//# filename:             TcpIpRoutines.h
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

#include "mbs_interface.h"
#include <string>

#ifndef TCPIPROUTINES_H
#define TCPIPROUTINES_H

enum SocketType {ST_server,ST_client};

class TCPIPSocket{

typedef _W64 unsigned int UINT_PTR, *PUINT_PTR;     //These definitions are made here locally by hand, even though they usually would be included via winsock2.h or windows.h.
typedef UINT_PTR        SOCKET;											//However, neither of these files can be included here because that would imply a preprocessor define
typedef unsigned char byte;													//#define ChooseColor ChooseColorA ... from the file CommDlg.h, which would mess with the function call
																										//mbs->ChooseColor in RigidBodyJoints.cpp, since the latter includes control.h, and with that, TcpIpRoutines.h in turn.

public:

	TCPIPSocket(){};

	TCPIPSocket(SocketType typei, char* ipdata=NULL, int errormsgflag=1, int infomsgflag=1);

	TCPIPSocket(SocketType typei, char* ipIn, unsigned short portIn, int errormsgflag=1, int infomsgflag=1);

	virtual void Set(SocketType typei, char* ipIn, unsigned short portIn, int errormsgflag=1, int infomsgflag=1);

	//initialize and configure sockets, establish connection; close prior connections if neccessary
	int SetUpConnection();

	//send a byte
	//return values cf. SendData(...)
	int SendByte(byte data);

	//receive a byte
	//return values cf. RecvData(...)
	int RecvByte(byte & data);

	//send a 4-byte integer number "data"
	//byte order converted to network byte order
	//return values cf. SendData(...)
	int SendInteger(int data);

	//receive a 4-byte integer number and store result in datain 
	//incoming byte order assumed to be network byte order
	//return values cf. RecvData(...)
	int RecvInteger(int & datain);

	//send "size" bytes data, starting from address "data" in the memory
	//return values: >0 ... number of successfully sent bytes
	//							 <0 ... an error occurred
	int SendData(char* data, int size);

	//receive "size" bytes data, starting at address "data" in the memory
	//return value: >0 ... number of successfully received bytes
	//					  	=0 ... s connection closed by server
	//							-1 ... an error occurred (e.g. connection time out, if a timeout was specified for the socket) 
	int RecvData(char* data, int size);

	//set a timeout for blocking receive-calls on socket TCPIPSocket::c; time is specified in ms;
	void SetRecvTimeOut(unsigned long time);

	//set a timeout for blocking send-calls on socket TCPIPSocket::c; time is specified in ms;
	void SetSendTimeOut(unsigned long time);

	//set a timeout for the blocking accept-call on the server socket TCPIPSocket::c; time is specified in ms;
	void SetAcceptTimeOut(unsigned long time);

	//set the per-socket buffer size for incoming data
	//void SetRecvBufferSize(int size);

	//set time delay for connect-loop
	void SetLoopTime(int time);

	//enable/disable output of error messages
	void SetErrorMessageMode(int flag=1);

	//enable/disable output of error messages
	void SetInfoMessageMode(int flag=1);

	//set TCPIPSocket::autoreconnect flag
	void SetAutoReconnect(int flag=0);

	//close connection & clean-up
	void CloseConnection();

	virtual ~TCPIPSocket();

protected:
	SocketType type;
	SOCKET s;  //server socket (->listen mode)
	SOCKET c;  //socket for data transfer with accepted client or client socket
	unsigned short port;
	char ipstring[16];
	int status; //contains status of socket s or c; 0 ... not initialized (ready for new connection), 1 ... connection open and working
	//0...OK, >0...an error occured, status is the number of the last error (WSAGetLastError())
	int autoreconnect; //1(0) ... (don't) attempt to establish a new connection with the given IP address and port automatically
										 //if an error occurred and previous connection was terminated; default is 0
	int looptime; //time delay for the connect-loop in ms; default is 1000 ms
	int errormsg; //1(0) flag to enable (disable) output of error messages; default is 0
	int infomsg;  //1(0) flag to enable (disable) output of info messages; default is 0
	//int recvbuffer; //size of buffer for receive calls in bytes, set via SetRecvBufferSize(...); if 0 - which is default - the system default for the buffer size is used
	int accepttimeout; //timeout for the blocking accept-call of a server-socket; disabled by default, and ignored for client-sockets; activated and set by SetAcceptTimeOut

	//Read TCP/IPv4 adress and port from file (path specified by "file"), or get manual user input if file not found
	void GetTCPIPConfig(char* file);

	//initialize winsock
	long StartWinsock();

	//create socket
	void CreateSocket();

	//bind server socket to given IP v4 adress and port
	long BindSocket();

	//set server socket into listen mode
	long ListenMode();

	//accept connection on server socket; TCPIPSocket::c then is the socket which is used for data transfer via that accepted connection
	int AcceptConnection();

	//start a loop in which, every TCPIPSocket::looptime milliseconds, a connection attempt is made by the client socket to the server socket
	long ConnectLoop();

	//as ConnectLoop, but just with one connection attempt
	long ConnectOnce();

	////manually "flush" buffer for incoming data
	//void FlushRecvBuffer();

	//everything which should be done on server/client side immediately after a connection has been established
	//to be overwritten in derived classes
	virtual void Initialize();

	//error handling; "err" and "errornumber" are error description and number
	//type: 0...fatal error --> program termination, 1...termination of TCP/IP connection (and possibly auto-reconnect)
	virtual int Error(int type, std::string err, int errornumber=0);

	//output of a message
	virtual void Message(std::string message);

	//output of an info message
	virtual void InfoMessage(std::string info);
};


class TCPIPHotInt: public TCPIPSocket{
public:

	TCPIPHotInt():TCPIPSocket(){mbs=NULL;};
	//TCPIPHotInt(SocketType typei, char* ipIn, unsigned short portIn, int errormsgflag=1, int infomsgflag=1):
	//TCPIPSocket(typei,ipIn,portIn,errormsgflag,infomsgflag){};

	virtual void Set(SocketType typei, char* ipIn, unsigned short portIn, int errormsgflag=1, int infomsgflag=1);
	void Set(MBS* mbsI, SocketType typei, char* ipIn, unsigned short portIn, int errormsgflag=1, int infomsgflag=1);

	virtual ~TCPIPHotInt(){};

protected:

	MBS* mbs;

	virtual int Error(int type, std::string err, int errornumber=0);
	virtual void Message(std::string message);
	virtual void InfoMessage(std::string info);

};

#endif