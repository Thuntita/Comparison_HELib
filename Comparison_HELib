#include <cstddef>
#include <sys/time.h>
#include "FHE.h"
#include "EncryptedArray.h"
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <gmp.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include "Ctxt.h"
#include "polyEval.h"
#include <algorithm>
#include <math.h>
#include "myUtils.h"

//define the numbers we want to compare at VECTOR_COUNT
#define VECTOR_COUNT 2
#define VECTOR_SIZE 1

//Equality recursive function (Z --> equal).
//Input the encrypted bits of a and b.
//Returns the encrypted boolean value of Enc(1) if a = b or Enc(0) otherwise.
Ctxt z(int i,int j, std::vector<Ctxt> digits1, std::vector<Ctxt> digits2, Ctxt enc1){
	// j is the length of integer bits.
	// i is the bits of integer.
	if (j == 1){
        	Ctxt computeZ1 = digits1[i];
        	computeZ1 += digits2[i];
        	computeZ1 += enc1;
       		return computeZ1;
		// computeZ1 return n-bits integer (2 = 0010 return 4 bits)
    	}else {
        	int l = ceil(j / 2);
       	 	Ctxt computeZ0 = z(i + l, j - l, digits1, digits2, enc1);
        	computeZ0 *= z(i, l, digits1, digits2, enc1);
        	return computeZ0;
		// computeZ0 return n-bits integer (2 = 0010 return 4 bits)                    
   	}
}

//Inequality recursive function (T --> greater than)
//Input the encrypted bits of a and b.
//Returns the encrypted boolean value of the comparison a > b.  
Ctxt t(int i,int j, std::vector<Ctxt> digits1, std::vector<Ctxt> digits2, Ctxt enc1){
    	if (j == 1){
        	Ctxt  computeT1 = digits1[i];
        	computeT1 *= digits2 [i];
    		computeT1 += digits1[i];
    		return computeT1;
    	}
    	else {
        	int l1 = ceil(j / 2);
        	Ctxt computeT0 = z(i + l1, j - l1, digits1, digits2, enc1); 
        	computeT0 *= t(i, l1, digits1, digits2, enc1); 
        	computeT0 += t(i + l1, j - l1, digits1, digits2, enc1);
        	return computeT0;                                    
    	}
}

//Minimum function with input two numbers (extracted to bits) with output the minimum of two numbers.
std::vector<Ctxt> selection(std::vector<Ctxt> digits2, std::vector<Ctxt> digits1, std::vector<Ctxt> digitsenc1, std::vector<Ctxt> digitst){
    for (int i = 0; i < digits1.size(); i++){
	// Return Minimum
        Ctxt reset = digitsenc1[0];
        digits1[i] *= digitst[0];
        digitsenc1[0] -= digitst[0];
        digits2[i] *= digitsenc1[0];
        digits1[i] += digits2[i];
        digitsenc1[0] = reset;
    }

    /*for (int i = 0; i < digits1.size(); i++){
	// Return Maximum
	Ctxt reset = digitsenc1[0];
	digitsenc1[0] += digitst[0];
	digits1[i] *= digitsenc1[0];
	digits2[i] *= digitst[0];
	digits1[i] += digits2[i];
	digitsenc1[0] = reset;
      }*/
    return digits1;
}

//function with input a list of encrypted numbers and output the minimum of them.
//Input an array list of encrypted integers.
//Output the encrypted minimum integer of the array.
std::vector<Ctxt> minimum(std::vector<Ctxt> list, Ctxt enc1){
    Ctxt ctxt1 = list[0];	// first plaintext

    // Use "extractDigits" function in Ctxt.h
    std::vector<Ctxt> digits1;
    extractDigits(digits1, ctxt1);

    for (int i = 1; i < list.size(); i++){
        std::vector<Ctxt> digitsenc1;
        extractDigits(digitsenc1, enc1);

        std::vector<Ctxt> digits2;
        extractDigits(digits2, list[i]);	// second plaintext
        
        Ctxt compute_t = t(0, digits1.size(), digits1, digits2, enc1);

        std::vector<Ctxt> digitst;
        extractDigits(digitst, compute_t);
        std::vector<Ctxt> min = selection(digits1, digits2, digitsenc1, digitst);
        digits1 = min;
    }
    return digits1;
}

int main(int argc, char **argv)
{
    /*** BEGIN INITIALIZATION ***/
    long m = 0;                   // Specific modulus
    long p = 2;                   // Plaintext base [default=2], should be a prime number
    long r = 4;                   // Lifting [default=1]
    long L = 16;                  // Number of levels in the modulus chain [default=heuristic]
    long c = 2;                   // Number of columns in key-switching matrix [default=2]
    long w = 5;                   // Hamming weight of secret key
    long d = 1;                   // Degree of the field extension [default=1]
    long k = 80;                  // Security parameter [default=80] 
    long s = 0;                   // Minimum number of slots [default=0]
    
    Timer tInit;
    tInit.start();
	
    std::cout << "Finding m... " << std::flush;
    m = FindM(k, L, c, p, d, s, 0);           // Find a value for m given the specified values
    
    std::cout << "m = " << m << std::endl;
	
    std::cout << "Initializing context... " << std::flush;
    FHEcontext context(m, p, r); 	      // Initialize context
    buildModChain(context, L, c);             // Modify the context, adding primes to the modulus chain
    std::cout << "OK!" << std::endl;

    std::cout << "Generating keys... " << std::flush;
    
    fstream secKeyFile("sk.txt", fstream::out|fstream::trunc);  
    assert(secKeyFile.is_open());
    writeContextBase(secKeyFile,context);
    secKeyFile << context << std::endl;

    FHESecKey sk(context);                    // Construct a secret key structure
    sk.GenSecKey(w);                          // Actually generate a secret key with Hamming weight
    addSome1DMatrices(sk);                    // Extra information for relinearization
    std::cout << "OK!" << std::endl;
    secKeyFile << sk << std::endl; 
    secKeyFile.close();

    /****INITIALIZATION END****/

    //Input numbers that we want to compare into vector 'u'
    long u[10]; 
    u[0] = 10;
    u[1] = 2;

    std::cout << "Plaintext: " << u[0] << std::endl; 
    std::cout << "Plaintext: " << u[1] << std::endl;

    //Encrypt numbers by using secret key.
    Timer tEnc;
    tEnc.start();
    Ctxt encU(sk), encV(sk), enc1(sk);
    sk.Encrypt(encU,to_ZZX(u[0]));
    sk.Encrypt(encV,to_ZZX(u[1]));
    sk.Encrypt(enc1,to_ZZX(1));
    tEnc.stop();
    std::cout << "Time for encryption: " << tEnc.elapsed_time() << "s" << std::endl;
    

    //Extract digits of encU and store them at vector:digitsU.    
    std::vector<Ctxt> digitsU;
    extractDigits(digitsU, encU);

    //Decrypt each digit of Ciphertext1 to see if extract digits work.
    long res[digitsU.size()];
    for (int i = 0; i < digitsU.size(); i++){
    	ZZX result;
    	sk.Decrypt(result, digitsU[i]);
	if (result[0] > (pow(p,r)) / 2){
		result[0] = result[0]-pow(p,r);
	}
    	res[i] = conv<long>(result[0]);
    }

    cout << "Binary1: ";
    size_t res_size = sizeof(res) / sizeof(res[0]);
    std::reverse(res, res + res_size);
    for (int i = 0; i < digitsU.size(); i++){
        std::cout << res[i];
    }
    std::cout << std::endl;

    //Extract digits of Ciphertext2 and store them at vector:digitsV.
    std::vector<Ctxt> digitsV;
    extractDigits(digitsV, encV);

    std::vector<Ctxt> digitsenc1;
    extractDigits(digitsenc1, enc1);

    //Decrypt each digit of encV to see if extractDigits work.
    long res1[digitsV.size()];
    for (int i = 0;i < digitsV.size(); i++){
    	ZZX result;
    	sk.Decrypt(result, digitsV[i]);
    	if (result[0] > (pow(p,r)) / 2){
		result[0] = result[0]-pow(p,r);
    	}
    	res1[i]=conv<long>(result[0]);
    }
    size_t res1_size = sizeof(res1) / sizeof(res1[0]);
    std::reverse(res1, res1 + res1_size);
    std::cout << "Binary2: ";
    for (int i = 0; i < digitsV.size(); i++){
        std::cout << res1[i];
    }
    std::cout << std::endl;

    //Input all encrypted numbers on list 
    std::vector<Ctxt> list;
    list.push_back(encU);
    list.push_back(encV);
    
    Timer timecompare;
    timecompare.start();
    //Calculate the encrypted minimum of the list.
    std::vector<Ctxt> Min = minimum(list, enc1);
    Ctxt ctxt1 = list[0];	// first plaintext

    // Use "extractDigits" function in Ctxt.h
    std::vector<Ctxt> digits1;
    extractDigits(digits1, ctxt1);

    Ctxt compute_t(sk);
    for (int i = 1; i < list.size(); i++){
        std::vector<Ctxt> digits2;
        extractDigits(digits2, list[i]);	// second plaintext
        
        compute_t = t(0, digits1.size(), digits1, digits2, enc1);

        std::vector<Ctxt> digitst;
        extractDigits(digitst, compute_t);
	
        std::vector<Ctxt> min = selection(digits1, digits2, digitsenc1, digitst);
        digits1 = min;
    }
    timecompare.stop();

    //Decryption of the minimum
    Timer tDec;
    tDec.start();
    long res2[digitsV.size()];
    for (int i = 0; i < digitsV.size(); i++){
        ZZX result2;
        sk.Decrypt(result2, digits1[i]);
        res2[i] = conv<long>(result2[0]);
    }
    tDec.stop();

    //Output the minimum number.
    std::cout<< "the minimum number is: ";
    int min = 0, flow;
    for (int i = 0; i < digitsV.size(); i++){
        flow = pow(2, i) * res2[i];
        min += flow;
    }
    std::cout << min << std::endl;
    std::cout << "the time to decryption is: " << tDec.elapsed_time() << "s" << std::endl;

    //The binary form of the number.
    std::cout << "Binary form is: " ;
    size_t res2_size = sizeof(res2)/sizeof(res2[0]);
    std::reverse(res2, res2 + res2_size);
    for (int i=0;i<digitsU.size();i++){
        std::cout << res2[i];
    }
    std::cout<< std::endl;
    std::cout << "the time to compute minimum is: " << timecompare.elapsed_time() << "s" << std::endl;
    std::cout << std::endl;
     return 0;
}
