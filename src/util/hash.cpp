#ifndef HASH_H
#define HASH_H
#include "hash.h"

Hash::Hash() {
    this->HASH_TABLE_SIZE = 30000000;
    this->table.resize(HASH_TABLE_SIZE, -1);
}

Hash::Hash(long table_size) {
    this->HASH_TABLE_SIZE = table_size;
    this->table.resize(this->HASH_TABLE_SIZE, -1);
}

unsigned int Hash::BKDRHash(char *key) {
    unsigned int seed = 131; // 31 131 1313 13131 131313 etc..
    unsigned int hash = 0;
    while (*key)
    {
        hash = hash * seed + (*key++);
    }
    return (hash % HASH_TABLE_SIZE);
}

void Hash::insert_key(char *key) {
    unsigned int pos = this->BKDRHash(key);
    while (this->table[pos] != -1)
        pos = (pos + 1) % this->HASH_TABLE_SIZE;
    this->table[pos] = this->keys.size();
    this->keys.push_back(strdup(key));
}

long Hash::search_key(char *key) {
    unsigned int pos = this->BKDRHash(key);
    while (1)
    {
        if (this->table[pos] == -1)
            return -1;
        if ( !strcmp(key, this->keys[ this->table[pos] ]) )
            return this->table[pos];
        pos = (pos + 1) % this->HASH_TABLE_SIZE;
    }
}

#endif
