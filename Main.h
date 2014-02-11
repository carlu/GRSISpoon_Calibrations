// Which parts of code to run:
#define SORT_CALIB 0
#define SORT_EFF   0
#define SORT_WAVES 0
#define SORT_PROP  1
#define SORT_DIFF  0

#define PRINT_OUTPUT 1
#define PRINT_FREQ 10000

#define MAX_EVENTS 0

#define DEBUG_TREE_LOOP 0


// ROOT Stuff
//#define ROOT_CACHE_SIZE 2000000000  //2GB
//#define ROOT_CACHE_SIZE 500000000  // 500MB
#define ROOT_VIRT_SIZE    500000000  //  500MB 
//                       /  /  /  /   

// Stuff about the experimental setup
#define CLOVERS  16
#define CRYSTALS  4
#define SEGS      8

#define EN_THRESH 2  // energies less than this keV ignored

#define PI 3.14159265359

// Structure to hold Mnemonic
struct Mnemonic	{
   int arrayposition;
   int	segment;
   std::string system;
   std::string subsystem;
   std::string arraysubposition;
   std::string collectedcharge;
   std::string outputsensor;
};

// Function to parse Mnemonic name:
void ParseMnemonic(std::string *name,Mnemonic *mnemonic);
int Col2Num(char Colour);
char Num2Col(int Crystal);


