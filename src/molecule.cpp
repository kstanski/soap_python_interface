#include "molecule.h"

#include <iostream>
#include <string.h>

Molecule **read_molecules(const char *filename, int molecules_no)
{
    /* create molecule array */
    Molecule **mol_arr = (Molecule **) malloc(molecules_no*sizeof(Molecule *));
    for (int i=0; i<molecules_no; i++)
    {
        mol_arr[i] = (Molecule *) malloc(sizeof(Molecule));
    }

    /* open file for reading */
    FILE *fp;
    fp = fopen(filename, "r");
    if(fp == NULL)
    {
        fprintf(stderr,"Error opening file\n");
        exit (1);
    }

    /* read line by line */
    char buf[256];
    int count_molecules = 0;
    int count_atoms = 0;
    Molecule *molecule = *mol_arr;   // points to the first Molecule.

    while (fgets(buf, sizeof(buf), fp) && count_molecules < molecules_no)
    {
        /* read words in a line */
        char *words[7]; // max 7 words in a line
        int word_count = 0;
        char *word;
        word = strtok(buf," ");
        while (word != NULL)
        {
            words[word_count] = word;
            word = strtok(NULL, " ");
            word_count++;
        }

        /* extract data from each line */
        switch (word_count)
        {
        case 1:
            count_atoms = 0;
            memset(molecule->types_total,0,ATOM_TYPES);
            molecule->atoms_no = atoi(words[0]);
            break;
        case 2:
            strncpy(molecule->id,words[0],ID_LEN);
            molecule->energy = atof(words[1]);
            break;
        case 7:
            molecule->atom_types[count_atoms] = type2index(words[0]);
            double ff_pos[DIMENSIONS], dft_pos[DIMENSIONS];
            for (int i=0; i<DIMENSIONS; i++)
            {
                ff_pos[i] = atof(words[i+1]);
                dft_pos[i] = atof(words[i+1+DIMENSIONS]);
            }
            molecule->ff_coords[count_atoms] = make_position(ff_pos);
            molecule->dft_coords[count_atoms] = make_position(dft_pos);
            molecule->types_total[type2index(words[0])]++;

            if (++count_atoms == molecule->atoms_no) // all atoms processed
            {
                count_molecules++;
                molecule = mol_arr[count_molecules]; // move pointer to the next Molecule.
            }
            break;
        default :
            fprintf(stderr,"Invalid entry\n");
        }
    }

    if (ferror(fp))
    {
        fprintf(stderr,"Error reading the file\n");
        exit (1);
    }
    return mol_arr;
}

Molecule **vectors2molecules(std::vector<intvector> atomic_no, std::vector<dmatrix> coords)
{
    /* create molecule array */
	int molecules_no = atomic_no.size();
    Molecule **mol_arr = (Molecule **) malloc(molecules_no*sizeof(Molecule *));
    for (int i=0; i<molecules_no; i++)
    {
        mol_arr[i] = (Molecule *) malloc(sizeof(Molecule));
    }

    int count_molecules = 0;

    while (count_molecules < molecules_no)
    {
        Molecule *molecule = mol_arr[count_molecules];
        for (int i=0; i<ATOM_TYPES; i++)
            molecule->types_total[i] = 0;

        intvector at_no = atomic_no[count_molecules];
        int atoms_no = at_no.size();
        molecule->atoms_no = atoms_no;
        int count_atoms = 0;
        while (count_atoms < atoms_no)
        {
            molecule->atom_types[count_atoms] = atomic_no2index(at_no[count_atoms]);
            double pos[DIMENSIONS];
            for (int i=0; i<DIMENSIONS; i++)
            {
                pos[i] = coords[count_molecules][count_atoms][i];
            }
            molecule->ff_coords[count_atoms] = make_position(pos);
            molecule->types_total[atomic_no2index(at_no[count_atoms])]++;
            count_atoms++;
        }

        count_molecules++;
    }

    return mol_arr;
}

Molecule *vectors2molecule(intvector atomic_no, dmatrix coords)
{
    /* create molecule array */
    Molecule *mol = (Molecule *) malloc(sizeof(Molecule));;

    for (int i=0; i<ATOM_TYPES; i++)
        mol->types_total[i] = 0;

    int atoms_no = atomic_no.size();
    mol->atoms_no = atoms_no;
    int count_atoms = 0;
    while (count_atoms < atoms_no)
    {
        mol->atom_types[count_atoms] = atomic_no2index(atomic_no[count_atoms]);
        double pos[DIMENSIONS];
        for (int i=0; i<DIMENSIONS; i++)
        {
            pos[i] = coords[count_atoms][i];
        }
        mol->ff_coords[count_atoms] = make_position(pos);
        mol->types_total[atomic_no2index(atomic_no[count_atoms])]++;
        count_atoms++;
    }

    return mol;
}

int type2index(char *type)
{
    if (strcmp("H",type) == 0) return 0;
    if (strcmp("C",type) == 0) return 1;
    if (strcmp("N",type) == 0) return 2;
    if (strcmp("O",type) == 0) return 3;
    if (strcmp("S",type) == 0) return 4;
    fprintf(stderr,"Unknown atom type\n");
    return -1;
}

int index2atomic_no(int idx)
{
    if (idx == 0) return 1;
    if (idx == 1) return 6;
    if (idx == 2) return 7;
    if (idx == 3) return 8;
    if (idx == 4) return 16;
    fprintf(stderr,"Unknown atom type\n");
    return -1;
}

int atomic_no2index(int at_no)
{
    if (at_no == 1) return 0;
    if (at_no == 6) return 1;
    if (at_no == 7) return 2;
    if (at_no == 8) return 3;
    if (at_no == 16) return 4;
    fprintf(stderr,"Unknown atom type\n");
    return -1;
}

Position make_position(double *val_arr)
{
    Position pos(val_arr[0],val_arr[1],val_arr[2]);
    return pos;
}

bool compare_molecules(Molecule *mol1, Molecule *mol2)
{
    return (*mol1).energy > (*mol2).energy;   //since energy is negative
}

int free_mol_array(Molecule **mol_arr, int molecules_no)
{
    for (int i=0; i<molecules_no; i++)
    {
        free(mol_arr[i]);
    }
    free(mol_arr);
    return 0;
}
