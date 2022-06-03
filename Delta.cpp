#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;


//create struct to save info for each probe
struct probe_info {
  int seq_index, pos;
  int fp_count, fn_count;
  double error_rate;
};

//create hash struct contanning 3 kind of hashes for each probe
struct Hash_info{
    long long hash;
    long long leftHash;
    long long rightHash;
};
//linked list to add to the table which each node
//contain probe_info struct
struct Node
{
    probe_info data;
    Node *next;
    Node (probe_info info, Node *n) { data = info; next = n; }
};

//declare basic variables to use in all functions
int N, N_delta, K = 100;
int prime = 10007;
int x = 7;
string *sequence;
bool *is_delta;

//create hash array for each probes
//which also contain left right and a whole hash
Hash_info **hashArray;

//create tables to save the hash values
//row is the number of sequences
//column equivalent to hash value
Node ***table;
Node ***leftTable;
Node ***rightTable;

//declare function prototype
probe_info eval_probe(int i, int p);
void read_input(void);
bool is_match(int i, int p, int j, int q,int size);
long long myhash(int i, int p,int size);
void loadHash( Hash_info **hashArray);
void allocateMemony();

int main(void)
{
    //read file to get N and N_delta
    read_input();
    
    //allocate memory for table, rightTable, and leftTable
    allocateMemony();
    
    //create hash array that contain hash of the whole probes,
    //half left and half right of probe
    hashArray = new Hash_info*[N];
    for(int a = 0; a < N; ++a)
    {
        hashArray[a] = new Hash_info[sequence[a].length() - K];
    }
    
    //calculate hash for all probes and add to the tables
    loadHash(hashArray);
     
    // Loop over all possible probes, remember the best one...
    probe_info best;
    best.error_rate = 99999;
    for (int i=0; i<N; i++) {
        cerr << ".";
        if (is_delta[i])
            for (int p=0; p<sequence[i].length()-K; p++)
            {
                // Evaluate probe at position p in sequence[i]
                probe_info info = eval_probe(i,p);
                if (info.error_rate < best.error_rate) best = info;
            }
    }
    
    //print out the output
    cerr << "\n";

    // Print out info about best
    cout << "this is the output for " << N << " sequences and " << "K = " << K << endl;
    cout << "Best probe: " << sequence[best.seq_index].substr(best.pos, K) << "\n";
    cout << "False positives: " << best.fp_count << "\n";
    cout << "False negatives: " << best.fn_count << "\n";
    cout << "Error_rate: " << best.error_rate << "\n";
    
    //free memory for tables
    for (int i = 0; i < N; ++i)
    {
        for(int j = 0; j < prime; ++j)
        {
            delete [] table[i][j];
            delete [] leftTable[i][j];
            delete [] rightTable[i][j];
        }
        delete [] table[i] ;
        delete [] leftTable[i];
        delete [] rightTable[i];
    }

    delete [] table ;
    delete [] leftTable;
    delete [] rightTable;
    
    //free memory for hashArray
    for(int a = 0; a < N; ++a)
    {
        delete [] hashArray[a];
    }
    delete [] hashArray;
    
    //free memory for sequences
    delete [] sequence;
    delete [] is_delta;
}
void allocateMemony()
{
    //allocate memory for table with
    //N row and column based on hash values
    table = new Node**[N];
    leftTable = new Node**[N];
    rightTable = new Node**[N];
    probe_info data;
    data.seq_index = 0;
    data.pos = 0;
    for (int i = 0; i < N; ++i)
    {
        table[i] = new Node*[prime];
        leftTable[i] = new Node*[prime];
        rightTable[i] = new Node*[prime];
        for(int j = 0; j < prime; ++j)
        {
            table[i][j] = new Node(data,NULL);
            leftTable[i][j] = new Node(data,NULL);
            rightTable[i][j] = new Node(data,NULL);
        }
    }
}

void read_input(void)
{
    ifstream input("covid.txt");
    string label, seq;
    while (input >> label >> seq) N++;
    input.clear();
    input.seekg(0);
    sequence = new string[N];
    is_delta = new bool[N];
    for (int i=0; i<N; i++)
    {
        input >> label >> sequence[i];
        is_delta[i] = label == "delta_variant";
        if (is_delta[i]) N_delta++;
    }
}

// Does position p in sequence[i] match position q in sequence[j] (1 mismatch char ok)?
bool is_match(int i, int p, int j, int q, int size)
{
    int mismatched = 0;
    for (int k=0; k< size; k++)
      if (sequence[i][p+k] != sequence[j][q+k]) {
        mismatched ++;
        //if more than 1 miss characters, return false
        if (mismatched > 1) return false;
      }
    //if 1 or less unmatched character, return true
    return true;
}

probe_info eval_probe(int i, int p)
{
    probe_info info;
    info.seq_index = i;
    info.pos = p;
    info.fp_count =  0;
    info.fn_count = 0;
    
    //compare the hash value of probe index i position j
    //too all the hash values saved in the table
    for (int j = 0; j < N; ++j)
    {
        
        if (is_delta[j]) info.fn_count ++;
        int count = 0;
        //find the matching probe based on hash of whole string
        for(Node *n = table[j][hashArray[i][p].hash]; n->next != NULL; n = n->next)
        {
            if( hashArray[i][p].leftHash == hashArray[j][n->data.pos].leftHash &&
               hashArray[i][p].rightHash == hashArray[j][n->data.pos].rightHash)
            {
                ++count;
                if (is_delta[j]) info.fn_count --;
                if (!is_delta[j]) info.fp_count ++;
                break;
            }
            
        }
        //if find a match, continue to next sequence
        if(count > 0) continue;
        
        //find the matching probe based on hash of left string
        for (Node *n = leftTable[j][hashArray[i][p].leftHash] ; n->next != NULL; n = n->next)
        {
            //if match hash value, compare the half right probes
            if(is_match(i, p + K/2, j, n->data.pos + K/2, round(K/2.0)))
            {
                ++count;
                if (is_delta[j]) info.fn_count --;
                if (!is_delta[j]) info.fp_count ++;
                //cout << "left" << endl;
                break;
            }
        }
        
        //if find a match, continue to next sequence
        if(count > 0) continue;
        
        //find the matching probe based on hash of right string
        for ( Node *n = rightTable[j][hashArray[i][p].rightHash]; n->next != NULL; n= n->next)
        {
            //if match hash value, compare the half left probes
            if(is_match(i, p, j, n->data.pos, K/2))
            {
                ++count;
                if (is_delta[j]) info.fn_count --;
                if (!is_delta[j]) info.fp_count ++;
                break;
            }
        }

    }
    //calculate error rate and save it
    double FPR = (double)info.fp_count / (N-N_delta);
    double FNR = (double)info.fn_count / N_delta;
    info.error_rate = 2.0 * FPR + 1.0 * FNR;
    return info;
}

void loadHash( Hash_info **hashArray)
{
    //calculate x^N-1 for hash
    long long pow1 = 1;
    for (int i = 0; i < K - 1 ; ++i)
    {
        pow1 = (pow1 * x)%prime;
    }
    //calculate x^N-1 for left hash
    long long pow2 = 1;
    for (int i = 0; i < K/2 - 1 ; ++i)
    {
        pow2 = (pow2 * x)%prime;
    }
    //calculate x^N-1 for right hash
    long long pow3 = 1;
    for (int i = 0; i < round(K/2.0) - 1 ; ++i)
    {
        pow3 = (pow3 * x)%prime;
    }
    
    for (int i = 0; i < N; ++i)
    {
        probe_info info;
        info.seq_index = i;
        info.pos = 0;
        
        //calculate first hash values each sequence
        long long hash = hashArray[i][0].leftHash = myhash(i,0,K/2);
        long long leftHash = hashArray[i][0].rightHash = myhash(i, K/2 ,round(K/2.0));
        long long rightHash = hashArray[i][0].hash = myhash(i, 0, K);
        //load those value to the table
        table[i][hashArray[i][0].hash] = new Node(info, table[i][hash]);
        leftTable[i][hashArray[i][0].leftHash] = new Node (info,leftTable[i][leftHash]);
        rightTable[i][hashArray[i][0].rightHash] = new Node (info , rightTable[i][rightHash]);
    
        for (int j = 1; j < sequence[i].length() - K; ++j)
        {
            //index and position for each probes
            info.seq_index = i;
            info.pos = j;
            
            //rolling hash for hash, left hash, and right hash
            hashArray[i][j].hash = (((hashArray[i][j-1].hash -
                                    ((long long)sequence[i][j-1]*pow1)%prime + prime) * x)% prime
                                    + sequence[i][j+K-1]) % prime;
            hashArray[i][j].leftHash = (((hashArray[i][j-1].leftHash -
                                       ((long long)sequence[i][j-1]*pow2)%prime + prime) * x)% prime
                                        + sequence[i][j+K/2-1]) % prime;
            hashArray[i][j].rightHash = (((hashArray[i][j-1].rightHash -
                                        ((long long)sequence[i][K/2+j-1]*pow3)%prime + prime) * x)% prime
                                        + sequence[i][j+K-1]) % prime;
            //get the equivalant hash value of each posiion
            long long hash = hashArray[i][j].hash;
            long long leftHash = hashArray[i][j].leftHash;
            long long rightHash = hashArray[i][j].rightHash;
            //load these value to the tables
            table[i][hash] = new Node(info, table[i][hash]);
            leftTable[i][leftHash] = new Node (info, leftTable[i][leftHash]);
            rightTable[i][rightHash] = new Node(info, rightTable[i][rightHash]);
            
        }
    }
}

long long myhash(int i, int p, int size)
{
    string s = sequence[i].substr(p,size);
    long long h = 0;
    //calculate hash and return it
    for (int i = 0; i < s.length();++i)
    {
        h = (((long long)h * x)%prime + s[i]) % prime;
    }
    return h;
}
