/* Exporting constants. See search.c for details */
/* TODO: Move overall design description here */
#ifndef SEARCH_H
#define SEARCH_H
/* Number of 'characters' that can fit in a long long. Calculated in search.c */
extern int chars_per_long;
#endif

/* 'safe overflow' approach to 'long' strings > chars_per_long in length */