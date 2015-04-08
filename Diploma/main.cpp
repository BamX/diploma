//
//  Copyright (c) 2015 BX23. All rights reserved.
//

#include <stdio.h>

#include "field.h"

int main(int argc, const char * argv[]) {
    Field *field = new Field(6, 3);
    //field->fill(75000.0);
    field->randomFill();
    field->print();
    field->transpose();
    field->print();

    return 0;
}
