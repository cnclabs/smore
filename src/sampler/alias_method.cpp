#include "alias_method.h"
#include <iostream>

AliasMethod::AliasMethod() {
}

AliasMethod::AliasMethod(std::vector<double>& distribution, const double power) {
    this->build(distribution, power);
}

void AliasMethod::build(std::vector<double>& distribution, const double power) {

    // normalization of vertices weights
    double sum, norm;
    std::vector<double> norm_prob;
    this->distribution_size = distribution.size();
    this->alias_position.resize(distribution.size());
    this->alias_probability.resize(distribution.size());
 
    sum = 0;
    std::vector<double>::iterator distribution_i;
    for (distribution_i=distribution.begin(); distribution_i!=distribution.end(); ++distribution_i)
    {
        sum += pow(*distribution_i, power);
    }
    norm = distribution.size()/sum;
    
    for (distribution_i=distribution.begin(); distribution_i!=distribution.end(); ++distribution_i)
    {
        norm_prob.push_back( pow(*distribution_i, power)*norm );
    }

    // block divison
    std::vector<long> small_block, large_block;
    
    for (long pos=0; pos!=norm_prob.size(); ++pos)
    {
        if ( norm_prob[pos]<1 )
        {
            small_block.push_back( pos );
        }
        else
        {
            large_block.push_back( pos );
        }
    }

    // assign alias table
    long small_pos, large_pos;

    while (small_block.size() && large_block.size())
    {
        small_pos = small_block.back();
        small_block.pop_back();
        large_pos = large_block.back();
        large_block.pop_back();

        this->alias_position[small_pos] = large_pos;
        this->alias_probability[small_pos] = norm_prob[small_pos];
        norm_prob[large_pos] = norm_prob[large_pos] + norm_prob[small_pos] - 1;
        if (norm_prob[large_pos] < 1)
        {
            small_block.push_back( large_pos );
        }
        else
        {
            large_block.push_back( large_pos );
        }
    }

    while (large_block.size())
    {
        large_pos = large_block.back();
        large_block.pop_back();
        this->alias_probability[large_pos] = 1.1; // any value > 1
    }

    while (small_block.size())
    {
        small_pos = small_block.back();
        small_block.pop_back();
        this->alias_probability[small_pos] = 1.1;
    }    
}

long AliasMethod::draw() {

    double sample_position = random_range(0, this->distribution_size);
    double sample_probability = random_range(0, 1);
    
    if (sample_probability < this->alias_probability[sample_position])
        return sample_position;
    else
        return this->alias_position[sample_position];
}

