//*******************************************************************************************************
//*******************************************************************************************************
//***********************************                          ******************************************
//***********************************        txHex2JSON        ******************************************
//***********************************                          ******************************************
//*******************************************************************************************************
//***************   -Converts a transaction hex to a JSON formatted data set-   *************************
//*******************************************************************************************************

#include <iostream>
#include <string>
#include <time.h>
#include <dirent.h>
#include <sstream>
#include <fstream>
#include <unistd.h>
#include "ensemble.h"
#include "PDBInterface.h"

vector<unsigned char> decodeTxHex(string hex);
vector<unsigned char> Parsehex(const char* hexchar);
static const signed char phexdigit[256] =
{ -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  0,1,2,3,4,5,6,7,8,9,-1,-1,-1,-1,-1,-1,
  -1,0xa,0xb,0xc,0xd,0xe,0xf,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,0xa,0xb,0xc,0xd,0xe,0xf,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, };


int main (int argc, char* argv[])
{
    //--Running parameters
    if (argc !=2)
    {
        cout << "txHex2JSON <transaction hex>" << endl;
        exit(1);
    }

    string hexString = argv[1];
    decodeTxHex(hexString);
    return 0;
}

vector<unsigned char> decodeTxHex(string hex)
{
    vector<unsigned char> txData = Parsehex(hex.c_str());

    for (UInt i = 0; i < txData.size(); i++)
    {
        stringstream ss;
        ss << txData[i] << hex;
        cout << ss << endl;
    }

    vector<unsigned char> reformattedtxData;// = TxToJSON(txData);

    return reformattedtxData;
}

vector<unsigned char> Parsehex(const char* hexchar)
{
    vector<unsigned char> vch;
    while (true)
    {
        while (isspace(*hexchar))
            hexchar++;
        signed char c = phexdigit[(unsigned char)*hexchar++];
        if (c == (signed char)-1)
            break;
        unsigned char n = (c << 4);
        c = phexdigit[(unsigned char)*hexchar++];
        if (c == (signed char)-1)
            break;
        n |= c;
        vch.push_back(n);
    }
    return vch;
}


/*vector<unsigned char> TxToJSON(vector<unsigned char> _txData)
{
    entry.push_back(Pair("txid", tx.GetHash().GetHex()));
    entry.push_back(Pair("version", tx.nVersion));
    entry.push_back(Pair("time", (boost::int64_t)tx.nTime));
    entry.push_back(Pair("locktime", (boost::int64_t)tx.nLockTime));
    Array vin;
    BOOST_FOREACH(const CTxIn& txin, tx.vin)
    {
        Object in;
        if (tx.IsCoinBase())
            in.push_back(Pair("coinbase", HexStr(txin.scriptSig.begin(), txin.scriptSig.end())));
        else
        {
            in.push_back(Pair("txid", txin.prevout.hash.GetHex()));
            in.push_back(Pair("vout", (boost::int64_t)txin.prevout.n));
            Object o;
            o.push_back(Pair("asm", txin.scriptSig.ToString()));
            o.push_back(Pair("hex", HexStr(txin.scriptSig.begin(), txin.scriptSig.end())));
            in.push_back(Pair("scriptSig", o));
        }
        in.push_back(Pair("sequence", (boost::int64_t)txin.nSequence));
        vin.push_back(in);
    }
    entry.push_back(Pair("vin", vin));
    Array vout;
    for (unsigned int i = 0; i < tx.vout.size(); i++)
    {
        const CTxOut& txout = tx.vout[i];
        Object out;
        out.push_back(Pair("value", ValueFromAmount(txout.nValue)));
        out.push_back(Pair("n", (boost::int64_t)i));
        Object o;
        ScriptPubKeyToJSON(txout.scriptPubKey, o, false);
        out.push_back(Pair("scriptPubKey", o));
        vout.push_back(out);
    }
    entry.push_back(Pair("vout", vout));
}*/

