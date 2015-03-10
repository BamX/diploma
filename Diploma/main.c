//
//  Copyright (c) 2015 BX23. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>

struct Field {
    size_t width, height;
    double **data;
};

struct Field *allocField(size_t width, size_t height) {
    struct Field *field = malloc(sizeof(*field));
    field->width = width;
    field->height = height;
    field->data = malloc(height * sizeof(double*));
    for (size_t i = 0; i < height; ++i) {
        field->data[i] = malloc(width * sizeof(double));
    }
    return field;
}

void freeField(struct Field *field) {
    for (size_t i = 0; i < field->height; ++i) {
        free(field->data[i]);
    }
    free(field->data);
    free(field);
}

void printField(struct Field *field) {
    for (size_t i = 0; i < field->height; ++i) {
        for (size_t j = 0; j < field->width; ++j) {
            printf("%.2f\t", field->data[i][j]);
        }
        printf("\n");
    }
}

void initFieldTemerature(struct Field *field, double temperature) {
    for (size_t i = 0; i < field->height; ++i) {
        for (size_t j = 0; j < field->width; ++j) {
            field->data[i][j] = temperature;
        }
    }
}

void passRow(struct Field *field, size_t row) {

}

void passColumn(struct Field *field, size_t col) {

}

int main(int argc, const char * argv[]) {
    struct Field *field = allocField(7, 7);
    initFieldTemerature(field, 75000.0);
    printField(field);
    freeField(field);

    return 0;
}
