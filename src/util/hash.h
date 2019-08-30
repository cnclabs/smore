#include <vector>
#include <string.h>

class Hash {
    /* Custom hash table
     */
    private:
        long HASH_TABLE_SIZE;
        unsigned int BKDRHash(char* key);

    public:
        Hash();
        Hash(long table_size);
        
        // variables
        std::vector< long > table; // init by -1
        std::vector< char* > keys;

        // operations
        void insert_key(char *key);
        long search_key(char *key);
};
