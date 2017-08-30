/*
ManyRig.h - Library facilitating the control of behaviral experiments.
Created by Duo Xu, September 24, 2016.
Released into the public domain.
*/

#ifndef ManyRig_h
#define ManyRig_h



#include "Arduino.h"

class ManyRig
{
public:
	ManyRig();

	void setDelimiter(char d);
	char getDelimiter();

	

	// Handling incoming messages
	void attachParser(void(*f)(void));
	void detachParser();
	unsigned int getInputIndex();
	unsigned long getInputValue();
	String getCmdString();

	void serialRead();
	void delay(unsigned long dur);
	bool delayUntil(bool(*f)(void));
	bool delayUntil(bool(*f)(void), unsigned long timeout);
	bool delayContinue(bool(*f)(void), unsigned long unitTime);



	// Sending formatted data message
	void attachSender(void(*f)(String));

	unsigned long sendData(char* tag, unsigned long t);
	unsigned long sendData(char* tag, unsigned long t, volatile byte num);
	unsigned long sendData(char* tag, unsigned long t, volatile int num);
	unsigned long sendData(char* tag, unsigned long t, volatile unsigned int num);
	unsigned long sendData(char* tag, unsigned long t, volatile long num);
	unsigned long sendData(char* tag, unsigned long t, volatile unsigned long num);
	unsigned long sendData(char* tag, unsigned long t, volatile float num);
	unsigned long sendData(char* tag, unsigned long t, volatile byte* dataArray, byte numData);
	unsigned long sendData(char* tag, unsigned long t, volatile int* dataArray, byte numData);
	unsigned long sendData(char* tag, unsigned long t, volatile unsigned int* dataArray, byte numData);
	unsigned long sendData(char* tag, unsigned long t, volatile long* dataArray, byte numData);
	unsigned long sendData(char* tag, unsigned long t, volatile unsigned long* dataArray, byte numData);
	unsigned long sendData(char* tag, unsigned long t, volatile float* dataArray, byte numData);



	// Control I/O
	unsigned long sendNumTTL(byte pin, unsigned long num);


protected:
	// Parsing
	char _delimiter = ',';
	unsigned int _numDelimiter = 0;
	unsigned long _inputVal = 0;
	String _cmdString = String();
	void (*_parserFunc)(void) = NULL;

	// Sending
	void (*_senderFunc)(String) = serialSend;

	static void serialSend(String);

};

#endif