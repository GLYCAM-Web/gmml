#include "../../includes/MolecularModeling/ring_shape_detection.hpp"

/*********************************************************************************************************/
/***************** Overview of port **********************************************************************/
// In May 2018 Oliver ported Spandana's C code for BFMP into GMML
// Method here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4210171/
// C program found by searching BFMP ring detection or probably broken by the time you read this:
// http://glycam.org/docs/othertoolsservice/download-docs/publication-materials/bfmp/
// I initially thought I would replace the classes with C++ and GMML classes , but after a few
// days I started just copying in the Glylib classes. The big change was to hard code canonicals.txt, as
// I don't think that will be changed very often, and how it was read in was awful. Hardcoding that
// meant I could delete a chunk of code. I also created a function does_conformation_match which
// somewhat reduces code replication, but boy this whole file is just tons of repeated code with
// different parameters. I ran out of time, and I'm thinking the more I change the harder it will
// be to find any bugs in this vs the original C code.

// There are a lot of geometry functions at the bottom. Some can be replaced by generic functions
// once they become available within gmml, others use the Glylib structs, so probably not unless
// you do a more extensive overhaul/replacment of the classes used here.
// Probably a namespace would be a good addition to hide the functions from other parts of GMML.
/***********************************************************************************************************/

constexpr auto PI = 3.14159265358979323846;


// This should not be in Assembly. Overloading to handle legacy code. Returning string so people can do as they wish.
std::string glylib::CalculateRingShapeBFMP( Glycan::Monosaccharide* mono )
{
  GeometryTopology::CoordinateVector ring_coordinates = gmml::GetCycleAtomCoordinates( mono );
  std::string bfmp = CalculateRingShapeBFMP(ring_coordinates);
  // OG is wondering why Monosaccharide isn't a class with Get'ers and Set'ers.
  mono->bfmp_ring_conformation_ = bfmp;
  return bfmp;
}

std::string glylib::CalculateRingShapeBFMP(GeometryTopology::CoordinateVector ring_coordinates, int cut_off)
{
    //****************************************************//
    //****************** canonicals.txt ******************//
    //****** OG is Hardcoding canonicals.txt here ********//

    int chair[12] = {1,2,4,5,0,2,3,5,0,1,3,4};
    int boat1[12] = {1,2,4,5,0,3,4,5,0,1,2,3};
    int boat2[12] = {0,1,3,4,2,3,4,5,0,1,2,5};
    int boat3[12] = {0,2,3,5,1,2,3,4,0,1,4,5};
    int skew1[8] = {1,2,3,5,0,2,4,5};
    int skew2[8] = {0,2,3,4,0,1,3,5};
    int skew3[8] = {1,3,4,5,0,1,2,4};
    int half1[4] = {2,3,4,5};
    int half2[4] = {0,3,4,5};
    int half3[4] = {0,1,4,5};
    int half4[4] = {0,1,2,5};
    int half5[4] = {0,1,2,3};
    int half6[4] = {1,2,3,4};

    //****************************************************//
    //****************** Output File *********************//
    fileset Oset;
    char *output = NULL;
    output=(char*)calloc(23,sizeof(char));
    sprintf(output,"ring_conformations.txt");
    Oset.N=strdup(output);
    Oset.F=myfopen(Oset.N,"w");
    //***************** End Output File ******************//

    std::stringstream result; // Will contain whatever normally goes into the output file.
    //result << "You: " << player << " CPU: " << cpu;
    //result.str();

    // std::cout << "Ok\n";
    double fifteen_dihedrals[16]; // hmmm.
    int list[60] = {0,1,2,3,0,1,2,4,0,1,2,5,0,1,3,4,0,1,3,5,0,1,4,5,0,2,3,4,0,2,3,5,0,2,4,5,0,3,4,5,1,2,3,4,1,2,3,5,1,2,4,5,1,3,4,5,2,3,4,5};
    int secondlist[30] = {4,5,3,5,3,4,2,5,2,4,2,3,1,5,1,4,1,3,1,2,0,5,0,4,0,3,0,2,0,1};
    char conformers[16] = {'q','t','q','d','t','q','t','d','t','q','q','t','d','t','q'};
    //StringVector conformers = {"q","t","q","d","t","q","t","d","t","q","q","t","d","t","q"};

    /*calculating the dihedral angles and planes*/
    double temp_dihedral0, temp_dihedral1, temp_dihedral2, temp_dihedral3;
    int j = 0, i = 0;
    for(i = 0; i < 15; i++)
    {
        temp_dihedral0 = calculateTorsionAngle(ring_coordinates.at(list[j]), ring_coordinates.at(list[j+1]), ring_coordinates.at(list[j+2]), ring_coordinates.at(list[j+3]));
        temp_dihedral1 = calculateTorsionAngle(ring_coordinates.at(list[j+1]), ring_coordinates.at(list[j+2]), ring_coordinates.at(list[j+3]), ring_coordinates.at(list[j]));
        temp_dihedral2 = calculateTorsionAngle(ring_coordinates.at(list[j+2]), ring_coordinates.at(list[j+3]), ring_coordinates.at(list[j]), ring_coordinates.at(list[j+1]));
        temp_dihedral3 = calculateTorsionAngle(ring_coordinates.at(list[j+3]), ring_coordinates.at(list[j]), ring_coordinates.at(list[j+1]), ring_coordinates.at(list[j+2]));
        //std::cout << temp_dihedral0 << ", " << temp_dihedral1 << ", " << temp_dihedral2 << ", " << temp_dihedral3 << "\n";
        fifteen_dihedrals[i] = (fabs(temp_dihedral0) + fabs(temp_dihedral1) + fabs(temp_dihedral2) + fabs(temp_dihedral3)) / 4 ;
        // printf("average_dihedrals are %lf\n",fifteen_dihedrals[i]);
        j = j + 4;
    }
    i = 0;
    j = 0;
    glylib::PlaneVector fifteen_planes;
    GeometryTopology::CoordinateVector four;
    double dihedral;
    int no = 4; // just no?
    for(i = 0; i < 15; i++)
    {
        //printf("Plane %d dihedral atoms are %d %d %d %d\n",i+1,list[j],list[j+1],list[j+2],list[j+3]);
        dihedral = calculateTorsionAngle(ring_coordinates.at(list[j]), ring_coordinates.at(list[j+1]), ring_coordinates.at(list[j+2]), ring_coordinates.at(list[j+3]));
        four.clear();
        four.push_back(ring_coordinates.at(list[j]));
        four.push_back(ring_coordinates.at(list[j+1]));
        four.push_back(ring_coordinates.at(list[j+2]));
        four.push_back(ring_coordinates.at(list[j+3]));

        fifteen_planes.push_back(get_plane_for_ring(no, four));
        //fifteen_planes[i]=get_plane_for_ring(no, four);
       // std::cout << "Plane: " << fifteen_planes.at(i).A << ", " << fifteen_planes.at(i).B << ", " << fifteen_planes.at(i).C << ", " << fifteen_planes.at(i).D << "\n";
        j = j + 4;
    }

    /*filtering all the planes with dihedral angles less than 10*/
    int no_planes=0;
    int bpsize=0;
    int sdsize=0;
    int *bestplanes;
    bestplanes=(int*)calloc(0,sizeof(int));
    int *sortedplanes;
    sortedplanes=(int*)calloc(0,sizeof(int));
    double *tendihedrals;
    tendihedrals=(double*)calloc(0,sizeof(double));
    double *sortdihedrals;
    sortdihedrals=(double*)calloc(0,sizeof(double));
    //double *threedihedrals; // commented out because not used
    //threedihedrals=(double*)calloc(3,sizeof(double));

    // cut_off = 10; // OG make this a passed variable with a 10 default
    int cut_off1 = (0 - cut_off);
    // printf("cutoffs are %d and %d\n", cut_off, cut_off1);

    for(i=0;i<15;i++){
        if(fifteen_dihedrals[i]<=cut_off && fifteen_dihedrals[i]>=cut_off1){
            bpsize++;
            bestplanes=(int*)realloc(bestplanes,(bpsize*sizeof(int)));
            bestplanes[bpsize-1]=i;
            tendihedrals=(double*)realloc(tendihedrals,(bpsize*sizeof(double)));
            tendihedrals[bpsize-1]=fifteen_dihedrals[i];
            sortdihedrals=(double*)realloc(sortdihedrals,(bpsize*sizeof(double)));
            sortdihedrals[bpsize-1]=fifteen_dihedrals[i];
            //   printf("%f\n",fabs(tendihedrals[bpsize-1]));
            no_planes++;
            //     printf("the planes are %d\n",i+1);
        }
    }

    double temp;
    for(i=0;i<bpsize;i++){
        for(j=i;j<bpsize;j++){
            if(fabs(sortdihedrals[i]) > fabs(sortdihedrals[j])){
                temp=sortdihedrals[i];
                sortdihedrals[i]=sortdihedrals[j];
                sortdihedrals[j]=temp;
            }
        }
    }

    for(i=0;i<bpsize;i++){
        for(j=0;j<bpsize;j++){
            if(sortdihedrals[i]==tendihedrals[j]){
                sdsize++;
                sortedplanes=(int*)realloc(sortedplanes,(sdsize*sizeof(int)));
                sortedplanes[sdsize-1]=bestplanes[j];
                //           printf("the sorted dihedrals are %f and planes are %d\n",sortdihedrals[i],bestplanes[j]);
            }
        }
    }
    //printf("sdsize is %d\n",sdsize);
    //    for(i=0;i<sdsize;i++){
    //        printf("sortedplanes are %d\n",sortedplanes[i]);
    //    }

    /*calculating the six dihedrals around the ring to check if it is a flat ring or an envelope*/
    j=0;
    //int six_list[25]= {1,2,3,4,2,3,4,5,3,4,5,6,4,5,6,1,5,6,1,2,6,1,2,3}
    int six_list[25]= {0,1,2,3,1,2,3,4,2,3,4,5,3,4,5,0,4,5,0,1,5,0,1,2};
    double *six_dihedrals;
    int sixcheck=0;
    six_dihedrals=(double*)calloc(7,sizeof(double));
    for(i=0;i<6;i++)
    {
        //dprint_coord_3D(&ring_coordinates.at(six_list[j]));
        //dprint_coord_3D(&ring_coordinates.at(six_list[j+1]));
        //dprint_coord_3D(&ring_coordinates.at(six_list[j+2]));
        //dprint_coord_3D(&ring_coordinates.at(six_list[j+3]));
        dihedral=calculateTorsionAngle(ring_coordinates.at(six_list[j]),ring_coordinates.at(six_list[j+1]),ring_coordinates.at(six_list[j+2]),ring_coordinates.at(six_list[j+3]));
        six_dihedrals[i] = dihedral * ( 180 / gmml::PI_RADIAN);
        if(six_dihedrals[i]<=5 && six_dihedrals[i]>=-5)
        {
            sixcheck++;
        }
        j=j+4;
    }

    if(sixcheck==6)
    {
        fprintf(Oset.F,"F\t");
        result << "F\t";
        //printf("This is a flat ring\n");
    }

    for(i=0;i<6;i++)
    {
        //printf("%lf\t",six_dihedrals[i]);
    }

    int envelope_check=0;
    int pre_e_check=0;
    //printf("%lf\t",fabs(fabs(six_dihedrals[0])-fabs(six_dihedrals[5])));
    if(fabs(fabs(six_dihedrals[0])-fabs(six_dihedrals[3])) <= 6.0){
        pre_e_check++;
    }
    //printf("%lf\t",fabs(fabs(six_dihedrals[1])-fabs(six_dihedrals[4])));
    if(fabs(fabs(six_dihedrals[1])-fabs(six_dihedrals[2])) <= 5.0){
        pre_e_check++;
    }
    //printf("%lf\t",fabs(fabs(six_dihedrals[2])-fabs(six_dihedrals[3])));
    if(fabs(fabs(six_dihedrals[4])-fabs(six_dihedrals[5])) <= 9.0){
        pre_e_check++;
    }
    //printf("pre check is %d\n",pre_e_check);

    if(pre_e_check<3){
        pre_e_check=0;
        if(fabs(fabs(six_dihedrals[1])-fabs(six_dihedrals[4])) <= 6.0){
            pre_e_check++;
        }
        if(fabs(fabs(six_dihedrals[2])-fabs(six_dihedrals[3])) <= 5.0){
            pre_e_check++;
        }
        if(fabs(fabs(six_dihedrals[5])-fabs(six_dihedrals[0])) <= 9.0){
            pre_e_check++;
        }
    }

    if(pre_e_check<3){
        pre_e_check=0;
        if(fabs(fabs(six_dihedrals[2])-fabs(six_dihedrals[5])) <= 6.0){
            pre_e_check++;
        }
        if(fabs(fabs(six_dihedrals[3])-fabs(six_dihedrals[4])) <= 5.0){
            pre_e_check++;
        }
        if(fabs(fabs(six_dihedrals[0])-fabs(six_dihedrals[1])) <= 9.0){
            pre_e_check++;
        }
    }


    if(pre_e_check<3){
        pre_e_check=0;
        if(fabs(fabs(six_dihedrals[0])-fabs(six_dihedrals[3])) <= 6.0){
            pre_e_check++;
        }
        if(fabs(fabs(six_dihedrals[4])-fabs(six_dihedrals[5])) <= 5.0){
            pre_e_check++;
        }
        if(fabs(fabs(six_dihedrals[1])-fabs(six_dihedrals[2])) <= 9.0){
            pre_e_check++;
        }
    }

    if(pre_e_check<3){
        pre_e_check=0;
        if(fabs(fabs(six_dihedrals[1])-fabs(six_dihedrals[4])) <= 6.0){
            pre_e_check++;
        }
        if(fabs(fabs(six_dihedrals[5])-fabs(six_dihedrals[0])) <= 5.0){
            pre_e_check++;
        }
        if(fabs(fabs(six_dihedrals[2])-fabs(six_dihedrals[3])) <= 9.0){
            pre_e_check++;
        }
    }

    if(pre_e_check<3){
        pre_e_check=0;
        if(fabs(fabs(six_dihedrals[2])-fabs(six_dihedrals[5])) <= 6.0){
            pre_e_check++;
        }
        if(fabs(fabs(six_dihedrals[0])-fabs(six_dihedrals[1])) <= 5.0){
            pre_e_check++;
        }
        if(fabs(fabs(six_dihedrals[3])-fabs(six_dihedrals[4])) <= 9.0){
            pre_e_check++;
        }
    }


    if(pre_e_check==3){
        envelope_check=1;
    }

    /*checking if there are two consecutive bad torsion angles in a row*/
    int *envelope_planes;
    envelope_planes=(int*)calloc(0,sizeof(int));
    int ssize;
    int esize=0;
    /*int *forty_planes;
   forty_planes=(int*)calloc(0,sizeof(int));
   int fosize=0;
   int esize=0;
   int ssize;
   int k=0;
   int fivecheck=0;
   if(sixcheck !=6 ){
           for(i=0;i<1;i++){
                   //printf("%lf and %lf\n",six_dihedrals[0],six_dihedrals[1]);
                   if((six_dihedrals[0] >=36.0 || six_dihedrals[0] <=-36.0) && (six_dihedrals[1]>=36.0 || six_dihedrals[1]<=-36.0)){
                           fivecheck++;
                           for(i=0;i<2;i++){
                                   fosize++;
                                   forty_planes=(int*)realloc(forty_planes,(fosize*sizeof(int)));
                                   forty_planes[fosize-1]=i;
                           }
                   }
                   //printf("%lf and %lf\n",six_dihedrals[1],six_dihedrals[2]);
                   if((six_dihedrals[1] >=36.0 || six_dihedrals[1] <=-36.0) && (six_dihedrals[2]>=36.0 || six_dihedrals[2]<=-36.0)){
                           fivecheck++;
                           for(i=1;i<3;i++){
                                   fosize++;
                                   forty_planes=(int*)realloc(forty_planes,(fosize*sizeof(int)));
                                   forty_planes[fosize-1]=i;
                           }

                   }
                   //printf("%lf and %lf\n",six_dihedrals[2],six_dihedrals[3]);
                   if((six_dihedrals[2] >=36.0 || six_dihedrals[2] <=-36.0) && (six_dihedrals[3]>=36.0 || six_dihedrals[3]<=-36.0)){
                           fivecheck++;
                           for(i=2;i<4;i++){
                                   fosize++;
                                   forty_planes=(int*)realloc(forty_planes,(fosize*sizeof(int)));
                                   forty_planes[fosize-1]=i;
                           }
                   }
                   //printf("%lf and %lf\n",six_dihedrals[3],six_dihedrals[4]);
                   if((six_dihedrals[3] >=36.0 || six_dihedrals[3] <=-36.0) && (six_dihedrals[4]>=36.0 || six_dihedrals[4]<=-36.0)){
                           fivecheck++;
                           for(i=3;i<5;i++){
                                   fosize++;
                                   forty_planes=(int*)realloc(forty_planes,(fosize*sizeof(int)));
                                   forty_planes[fosize-1]=i;
                           }

                   }
                   //printf("%lf and %lf\n",six_dihedrals[4],six_dihedrals[5]);
                   if((six_dihedrals[4] >=36.0 || six_dihedrals[4] <=-36.0) && (six_dihedrals[5]>=36.0 || six_dihedrals[5]<=-36.0)){
                           fivecheck++;
                           for(i=4;i<6;i++){
                                   fosize++;
                                   forty_planes=(int*)realloc(forty_planes,(fosize*sizeof(int)));
                                   forty_planes[fosize-1]=i;
                           }
                   }
                   //printf("%lf and %lf\n",six_dihedrals[5],six_dihedrals[0]);
                   if((six_dihedrals[5] >=36.0 || six_dihedrals[5] <=-36.0) && (six_dihedrals[0]>=36.0 || six_dihedrals[0]<=-36.0)){
                           fivecheck++;
                           fosize++;
                           forty_planes=(int*)realloc(forty_planes,(fosize*sizeof(int)));
                           forty_planes[fosize-1]=5;
                           fosize++;
                           forty_planes=(int*)realloc(forty_planes,(fosize*sizeof(int)));
                           forty_planes[fosize-1]=0;

                   }
           }
   }
//printf("fivecheck is %d\n",fivecheck);

if(fivecheck==1){
   //printf("This is an envelope\n");
}

for(i=0;i<fosize;i++){
   //printf("The dihedrals are %d\n",forty_planes[i]);
}

int *thirty_planes;
thirty_planes=(int*)calloc(0,sizeof(int));
int tsize=0;
int focheck=0;
if(fivecheck==1){
   for(i=0;i<6;i++){
           //printf("forty_planes[0] is %d and forty_planes[1] is %d\n",forty_planes[0],forty_planes[1]);
           if(i!= forty_planes[0] && i!= forty_planes[1]){
                   //printf("inside if %d\n",i);
                   if(six_dihedrals[i]>30.0 || six_dihedrals[i]<-30.0){
                           tsize++;
                           thirty_planes=(int*)realloc(thirty_planes,(tsize*sizeof(int)));
                           thirty_planes[tsize-1]=i;
                           focheck++;

                   }
           }
   }
}

//printf("thirty planes are %d\n",thirty_planes[0]);

//printf("focheck is %d\n",focheck);
int thcheck=0;

if(focheck ==1){
   //printf("inside focheck\n");
   for(i=0;i<6;i++){
           if(i != forty_planes[0] && i!= forty_planes[1]){
                   if(i!=thirty_planes[0]){
                           //printf("inside if\n");
                           if((six_dihedrals[i]>=20.0 && six_dihedrals[i]<=30.0) || (six_dihedrals[i]<=-20.0 && six_dihedrals[i]>=-30.0)){
                                   thcheck=1;
                           }
                   }
           }
   }
}

//printf("thcheck is %d\n",thcheck);*/
    /*checking if the opposite angles are less than zero*/
    int fiveplanes=0;
    int fivecheck2=0;
    if(envelope_check==1){
        for(i=0;i<1;i++){
            //printf("%lf and %lf\n",six_dihedrals[0],six_dihedrals[1]);
            if((six_dihedrals[0] <=12.0 && six_dihedrals[0] >=-12.0) && (six_dihedrals[1]<=12.0 && six_dihedrals[1]>=-12.0)){
                fivecheck2=1;
                ssize=0;
                for(j=ssize;j<ssize+8;j++){
                    esize++;
                    envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                    envelope_planes[esize-1]=six_list[j];

                }
                //printf("yes\n");
                fiveplanes++;
            }
            //printf("%lf and %lf\n",six_dihedrals[1],six_dihedrals[2]);
            if((six_dihedrals[1] <=12.0 && six_dihedrals[1] >=-12.0) && (six_dihedrals[2]<=12.0 && six_dihedrals[2]>=-12.0)){
                fivecheck2=1;
                ssize=4;
                for(j=ssize;j<ssize+8;j++){
                    esize++;
                    envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                    envelope_planes[esize-1]=six_list[j];

                }
                //printf("yes\n");
                fiveplanes++;
            }
            //printf("%lf and %lf\n",six_dihedrals[2],six_dihedrals[3]);
            if((six_dihedrals[2] <=12.0 && six_dihedrals[2] >=-12.0) && (six_dihedrals[3]<=12.0 && six_dihedrals[3]>=-12.0)){
                fivecheck2=1;
                ssize=8;
                for(j=ssize;j<ssize+8;j++){
                    esize++;
                    envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                    envelope_planes[esize-1]=six_list[j];

                }
                //printf("yes\n");
                fiveplanes++;
            }
            //printf("%lf and %lf\n",six_dihedrals[3],six_dihedrals[4]);
            if((six_dihedrals[3] <=12.0 && six_dihedrals[3] >=-12.0) && (six_dihedrals[4]<=12.0 && six_dihedrals[4]>=-12.0)){
                fivecheck2=1;
                ssize=12;
                for(j=ssize;j<ssize+8;j++){
                    esize++;
                    envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                    envelope_planes[esize-1]=six_list[j];

                }
                //printf("yes\n");
                fiveplanes++;
            }
            //printf("%lf and %lf\n",six_dihedrals[4],six_dihedrals[5]);
            if((six_dihedrals[4] <=12.0 && six_dihedrals[4] >=-12.0) && (six_dihedrals[5]<=12.0 && six_dihedrals[5]>=-12.0)){
                fivecheck2=1;
                ssize=16;
                for(j=ssize;j<ssize+8;j++){
                    esize++;
                    envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                    envelope_planes[esize-1]=six_list[j];

                }
                //printf("yes\n");
                fiveplanes++;
            }
            //printf("%lf and %lf\n",six_dihedrals[5],six_dihedrals[0]);
            if((six_dihedrals[5] <=12.0 && six_dihedrals[5] >=-12.0) && (six_dihedrals[0]<=12.0 && six_dihedrals[0]>=-12.0)){
                fivecheck2=1;
                ssize=20;
                for(j=ssize;j<ssize+4;j++){
                    esize++;
                    envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                    envelope_planes[esize-1]=six_list[j];

                }
                ssize=0;
                for(j=ssize;j<ssize+4;j++){
                    esize++;
                    envelope_planes=(int*)realloc(envelope_planes,(esize*sizeof(int)));
                    envelope_planes[esize-1]=six_list[j];

                }
                //printf("yes\n");
                fiveplanes++;
            }
        }

    }


    //printf("five planes is %d\n",fiveplanes*2);

    /*for(i=0;i<esize;i++){
   printf("%d\n",envelope_planes[i]);
}*/


    /*determining the conformation of the envelope*/
    double *d;
    d=(double*)calloc(6,sizeof(double));
    i=0;
    j=0;
    int k = 0;
    int echeck=0;
    int fsize=0;
    int *fiveatom_planes;
    fiveatom_planes=(int*)calloc(0,sizeof(int));

    if(envelope_check==1 && fivecheck2==1){
        //printf("fiveatomplane\n");
        while(i<fiveplanes*2*4){
            j=i;
            while(j<i+4){
                echeck=0;
                for(k=i+4;k<i+8;k++){
                    if(envelope_planes[j] == envelope_planes[k]){
                        //printf("%d\n",envelope_planes[k]);
                        echeck=1;
                    }
                }
                if(echeck==0){
                    fsize++;
                    fiveatom_planes=(int*)realloc(fiveatom_planes,(fsize*sizeof(fiveatom_planes)));
                    fiveatom_planes[fsize-1]=envelope_planes[j];
                    for(k=i+4;k<i+8;k++){
                        fsize++;
                        fiveatom_planes=(int*)realloc(fiveatom_planes,(fsize*sizeof(fiveatom_planes)));
                        fiveatom_planes[fsize-1]=envelope_planes[k];
                    }
                    j=i+4;
                    //printf("%d\n",envelope_planes[j]);
                }

                j++;
            }
            i=i+8;
        }//while(i<fiveplanes*2*4)


        //        for(i=0;i<fsize;i++){
        //            printf("%d\n",fiveatom_planes[i]);
        //        }

        int ringatoms[6]= {0,1,2,3,4,5};
        ssize=0;
        int *sixth_atom;
        sixth_atom=(int*)calloc(0,sizeof(sixth_atom));
        i=0;
        while(i<fiveplanes*5){
            j=i;
            k=0;
            while(k<6){
                echeck=0;
                for(j=i;j<i+5;j++){
                    if(fiveatom_planes[j]==ringatoms[k]){
                        echeck=1;
                    }
                }

                if(echeck==0){
                    //printf("%d\n",ringatoms[k]);
                    ssize++;
                    sixth_atom=(int*)realloc(sixth_atom,(ssize*sizeof*(sixth_atom)));
                    sixth_atom[ssize-1]=ringatoms[k];

                }
                k++;
            }

            i=i+5;
        }//while (i<fiveplanes*5)

        GeometryTopology::CoordinateVector five;
        //coord_3D **five;
        //five=(coord_3D**)calloc(6,sizeof(coord_3D*));
        plane *fiveatomplane;
        fiveatomplane=(plane*)calloc(1,sizeof(plane));
        int n = 5;



        for(i=0;i<fiveplanes;i++){
            five.clear();
            five.push_back(ring_coordinates.at(fiveatom_planes[i]));
            five.push_back(ring_coordinates.at(fiveatom_planes[i+1]));
            five.push_back(ring_coordinates.at(fiveatom_planes[i+2]));
            five.push_back(ring_coordinates.at(fiveatom_planes[i+3]));
            five.push_back(ring_coordinates.at(fiveatom_planes[i+4]));
            fiveatomplane[0]=get_plane_for_ring(n,five);
            d[i]=get_signed_distance_from_point_to_plane(fiveatomplane[0],ring_coordinates.at(sixth_atom[i]));
            //printf("distances are %f\n",d[i]);

        }

        double min=0.0;
        int min_atom=0;
        min=d[0];
        min_atom=sixth_atom[0];
        for(i=0;i<fiveplanes;i++)
        {
            if(d[i]<min)
            {
                min=d[i];
                min_atom=sixth_atom[i];
            }
        }
        if(min<0.0)
        {
            fprintf(Oset.F,"\tE%d\t", min_atom+1);
            fprintf(Oset.F,"\tp%d\t", min_atom+1);
            result << "\tE" << min_atom+1 << "\t\tp" << min_atom+1 << "\t";
        }
        else
        {
            fprintf(Oset.F,"\t%dE\t", min_atom+1);
            fprintf(Oset.F,"\t%dp\t", min_atom+1);
            result << "\t" << min_atom+1 << "E\t\t" << min_atom+1 << "p\t";
        }

        //	printf("The min distance is %f and the atom is %d\n",min,min_atom+1);

    }

    //printf("%d and %d \n",fivecheck,fivecheck2);
    if(fivecheck2==0)
    {
        //printf("The number of good planes are %d\n",no_planes);
    }
    if(fivecheck2==0 && no_planes==0)
    {
        fprintf(Oset.F,"\t-\t\t\tm\t");
        result << "\t-\t\t\tm\t";
    }
    int ccheck=0, bcheck=0;
    if(fivecheck2==0 && no_planes>0)
    {
        if(sdsize>=3)
        {
            if(does_conformation_match(3, sortedplanes, chair, 12, list, 3))
            {
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[12],ring_coordinates.at(3));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[12],ring_coordinates.at(0));
                //printf("the distances are %f and %f\n",d[0],d[1]);
                if(d[0]<0.0 && d[1] > 0.0)
                {
                    fprintf(Oset.F,"\t1C4\t");
                    result << "\t1C4\t";
                }
                if(d[0]>0.0 && d[1] <0.0)
                {
                    fprintf(Oset.F,"\t4C1\t");
                    result << "\t4C1\t";
                }
                ccheck=1;
            }
        }
        //printf("ccheck is %d\n",ccheck);
        if(ccheck==0 && sdsize >=3){
            if(does_conformation_match(sdsize, sortedplanes, boat1, 12, list, 3))
            {
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[12],ring_coordinates.at(3));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[12],ring_coordinates.at(0));
                if(d[0]>0.0 && d[1]>0.0)
                {
                    bcheck=1;
                    fprintf(Oset.F,"\t14B\t");
                    result << "\t14B\t";
                }
                if(d[0]<0.0 && d[1]<0.0)
                {
                    bcheck=1;
                    fprintf(Oset.F,"\tB14\t");
                    result << "\tB14\t";
                }
            }
            if(does_conformation_match(sdsize, sortedplanes, boat2, 12, list, 3))
            {
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[3],ring_coordinates.at(2));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[3],ring_coordinates.at(5));
                if(d[0]>0.0 && d[1]>0.0)
                {
                    fprintf(Oset.F,"\tO3B\t");
                    result << "\tO3B\t";
                    bcheck=1;
                }
                if(d[0]<0.0 && d[1]<0.0)
                {
                    fprintf(Oset.F,"\tBO3\t");
                    result << "\tBO3\t";
                    bcheck=1;
                }
            }
            if(does_conformation_match(sdsize, sortedplanes, boat3, 12, list, 3))
            {
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[7],ring_coordinates.at(1));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[7],ring_coordinates.at(4));
                if(d[0]>0.0 && d[1]>0.0)
                {
                    fprintf(Oset.F,"\t25B\t");
                    result << "\t25B\t";
                    bcheck=1;
                }
                if(d[0]<0.0 && d[1]<0.0)
                {
                    fprintf(Oset.F,"\tB25\t");
                    result << "\tB25\t";
                    bcheck=1;
                }
            }
        }
        //printf("ccheck and bcheck are %d and %d\n",bcheck,ccheck);
        int scheck = 0;
        int nsize = sdsize;
        if(ccheck==0 && bcheck==0)
        {
            if(nsize >= 3)
            {
                nsize = 3;
            }
            if(does_conformation_match(nsize, sortedplanes, skew1, 8, list, 2))
            {
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[11],ring_coordinates.at(0));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[11],ring_coordinates.at(4));
                if(d[0]<0.0 && d[1]>0.0)
                {
                    scheck=1;
                    fprintf(Oset.F,"\t5S1\t");
                    result << "\t5S1\t";
                }
                if(d[0]>0.0 && d[1]<0.0)
                {
                    scheck=1;
                    fprintf(Oset.F,"\t1S5\t");
                    result << "\t1S5\t";
                }
            }
            if(does_conformation_match(nsize, sortedplanes, skew2, 8, list, 2))
            {
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[6],ring_coordinates.at(5));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[6],ring_coordinates.at(1));
                if(d[0]<0.0 && d[1]>0.0)
                {
                    scheck=1;
                    fprintf(Oset.F,"\t2SO\t");
                    result << "\t2SO\t";
                }
                if(d[0]>0.0 && d[1]<0.0)
                {
                    fprintf(Oset.F,"\tOS2\t");
                    result << "\tOS2\t";
                }
            }
            if(does_conformation_match(nsize, sortedplanes, skew3, 8, list, 2))
            {
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[13],ring_coordinates.at(2));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[13],ring_coordinates.at(0));
                if(d[0]<0.0 && d[1]>0.0)
                {
                    scheck=1;
                    fprintf(Oset.F,"\t1S3\t");
                    result << "\t1S3\t";
                }
                if(d[0]>0.0 && d[1]<0.0)
                {
                    scheck=1;
                    fprintf(Oset.F,"\t3S1\t");
                    result << "\t3S1\t";
                }
            }
        }
        //printf("bcheck,ccheck and scheck are %d,%d and %d\n",bcheck,ccheck,scheck);
        int hcheck=0;
        if(ccheck==0 && bcheck==0 && scheck==0)
        {
            if(does_conformation_match(1, sortedplanes, half1, 4, list, 1))
            {
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[14],ring_coordinates.at(0));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[14],ring_coordinates.at(1));
                if(d[0]<0.0 && d[1]>0.0)
                {
                    //OGi=Cslurp.n;
                    hcheck=1;
                    fprintf(Oset.F,"\t2H1\t");
                    result << "\t2H1\t";
                }
                if(d[0]>0.0 && d[1]<0.0)
                {
                    //OGi=Cslurp.n;
                    hcheck=1;
                    fprintf(Oset.F,"\t1H2\t");
                    result << "\t1H2\t";
                }
            }
            if(does_conformation_match(1, sortedplanes, half2, 4, list, 1))
            {
                //printf("inside Half2\n");
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[9],ring_coordinates.at(2));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[9],ring_coordinates.at(1));
                //printf("%f and %f\n",d[0],d[1]);
                if(d[0]<0.0 && d[1]>0.0)
                {
                    //OGi=Cslurp.n;
                    hcheck=1;
                    fprintf(Oset.F,"\t2H3\t");
                    result << "\t2H3\t";
                }
                if(d[0]>0.0 && d[1]<0.0)
                {
                    //OGi=Cslurp.n;
                    hcheck=1;
                    fprintf(Oset.F,"\t3H2\t");
                    result << "\t3H2\t";
                }
            }
            if(does_conformation_match(1, sortedplanes, half3, 4, list, 1))
            {
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[5],ring_coordinates.at(2));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[5],ring_coordinates.at(3));
                if(d[0]<0.0 && d[1]>0.0)
                {
                    //OGi=Cslurp.n;
                    hcheck=1;
                    fprintf(Oset.F,"\t4H3\t");
                    result << "\t4H3\t";
                }
                if(d[0]>0.0 && d[1]<0.0)
                {
                    //OGi=Cslurp.n;
                    hcheck=1;
                    fprintf(Oset.F,"\t3H4\t");
                    result << "\t3H4\t";
                }
            }
            if(does_conformation_match(1, sortedplanes, half4, 4, list, 1))
            {
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[2],ring_coordinates.at(4));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[2],ring_coordinates.at(3));
                if(d[0]<0.0 && d[1]>0.0)
                {
                    //OGi=Cslurp.n;
                    hcheck=1;
                    fprintf(Oset.F,"\t4H5\t");
                    result << "\t4H5\t";
                }
                if(d[0]>0.0 && d[1]<0.0)
                {
                    //OGi=Cslurp.n;
                    hcheck=1;
                    fprintf(Oset.F,"\t5H4\t");
                    result << "\t5H4\t";
                }
            }
            if(does_conformation_match(1, sortedplanes, half5, 4, list, 1))
            {
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[0],ring_coordinates.at(4));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[0],ring_coordinates.at(5));
                if(d[0]<0.0 && d[1]>0.0)
                {
                    //OGi=Cslurp.n;
                    hcheck=1;
                    fprintf(Oset.F,"\tOh5\t");
                    result << "\tOh5\t";
                }
                if(d[0]>0.0 && d[1]<0.0)
                {
                    //OGi=Cslurp.n;
                    hcheck=1;
                    fprintf(Oset.F,"\t5HO\t");
                    result << "\t5HO\t";
                }
            }
            if(does_conformation_match(1, sortedplanes, half6, 4, list, 1))
            {
                d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[10],ring_coordinates.at(0));
                d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[10],ring_coordinates.at(5));
                if(d[0]<0.0 && d[1]>0.0)
                {
                    //OGi=Cslurp.n;
                    hcheck=1;
                    fprintf(Oset.F,"\tOH1\t");
                    result << "\tOH1\t";
                }
                if(d[0]>0.0 && d[1]<0.0)
                {
                    ////OGi=Cslurp.n;
                    hcheck=1;
                    fprintf(Oset.F,"\t1HO\t");
                    result << "\t1HO\t";
                }
            }
        }
        if(ccheck==0 && bcheck==0 && scheck==0 && hcheck==0)
        {
          //  printf("Bingo\n");
            fprintf(Oset.F,"\t-\t");
            result << "\t-\t";
        }
        for(i=0;i<no_planes;i++)
        {
            //printf("sortedplanes are %d\n",sortedplanes[i]);
            //printf("conformers are %s\n",conformers[sortedplanes[i]]);
            //printf("secondlistatoms are %d && %d\n",secondlist[sortedplanes[i]*2],secondlist[sortedplanes[i]*2+1]);
            d[0]=get_signed_distance_from_point_to_plane(fifteen_planes[sortedplanes[i]],ring_coordinates.at(secondlist[sortedplanes[i]*2]));
            d[1]=get_signed_distance_from_point_to_plane(fifteen_planes[sortedplanes[i]],ring_coordinates.at(secondlist[sortedplanes[i]*2+1]));
            //printf("%f and %f\n",d[0],d[1]);
            if(d[0]<0.0 && d[1]>0.0)
            {
                fprintf(Oset.F,"\t%d%c%d(%f)\t", secondlist[sortedplanes[i]*2+1]+1, conformers[sortedplanes[i]], secondlist[sortedplanes[i]*2]+1, fifteen_dihedrals[sortedplanes[i]]);
                result << "\t" << secondlist[sortedplanes[i]*2+1]+1 << conformers[sortedplanes[i]] << secondlist[sortedplanes[i]*2]+1 << "(" << fifteen_dihedrals[sortedplanes[i]] << ")" << "\t";
            }
            if(d[0]>0.0 && d[1]<0.0)
            {
                fprintf(Oset.F,"\t%d%c%d(%f)\t", secondlist[sortedplanes[i]*2]+1, conformers[sortedplanes[i]], secondlist[sortedplanes[i]*2+1]+1, fifteen_dihedrals[sortedplanes[i]]);
                result << "\t" << secondlist[sortedplanes[i]*2]+1 << conformers[sortedplanes[i]] << secondlist[sortedplanes[i]*2+1]+1 << "(" << fifteen_dihedrals[sortedplanes[i]] << ")" << "\t";
            }
            if(d[0]>0.0 && d[1] >0.0)
            {
                fprintf(Oset.F,"\t%d%d%c(%f)\t", secondlist[sortedplanes[i]*2]+1, secondlist[sortedplanes[i]*2+1]+1, conformers[sortedplanes[i]], fifteen_dihedrals[sortedplanes[i]]);
                result << "\t" << secondlist[sortedplanes[i]*2]+1 << secondlist[sortedplanes[i]*2+1]+1 << conformers[sortedplanes[i]] << "(" << fifteen_dihedrals[sortedplanes[i]] << ")" << "\t";
            }
            if(d[0]<0.0 && d[1]<0.0)
            {
                fprintf(Oset.F,"\t%c%d%d(%f)\t",conformers[sortedplanes[i]],secondlist[sortedplanes[i]*2]+1,secondlist[sortedplanes[i]*2+1]+1,fifteen_dihedrals[sortedplanes[i]]);
                result << "\t" << conformers[sortedplanes[i]] << secondlist[sortedplanes[i]*2]+1 << secondlist[sortedplanes[i]*2+1]+1 << "(" << fifteen_dihedrals[sortedplanes[i]] << ")" <<  "\t";
            }
        }
    } //fivecheck2==0
    fprintf(Oset.F,"\n");
   // printf("\n");
   // printf("OUTPUT FILE:ring_conformations.txt\n");
    result << "\n";
    return result.str();
}

// Copied this from gmml Assembly. Should use a generic one from Geometry once it becomes available
double glylib::calculateTorsionAngle(GeometryTopology::Coordinate *coord1, GeometryTopology::Coordinate *coord2, GeometryTopology::Coordinate *coord3, GeometryTopology::Coordinate *coord4)
{
    double current_dihedral = 0.0;

    // Oliver updates to solve memory leaks
    GeometryTopology::Coordinate b1 = *coord2; // deep copy
    GeometryTopology::Coordinate b2 = *coord3;
    GeometryTopology::Coordinate b3 = *coord4;
    GeometryTopology::Coordinate b4 = b2;
    b1.operator -(*coord1);
    b2.operator -(*coord2);
    b3.operator -(*coord3);
    b4.operator *(-1);

    GeometryTopology::Coordinate b2xb3 = b2; // deep copy
    b2xb3.CrossProduct(b3);

    GeometryTopology::Coordinate b1_m_b2n = b1;
    b1_m_b2n.operator *(b2.length());

    GeometryTopology::Coordinate b1xb2 = b1; // deep copy
    b1xb2.CrossProduct(b2);

    current_dihedral = atan2(b1_m_b2n.DotProduct(b2xb3), b1xb2.DotProduct(b2xb3));
    return (current_dihedral * (180 / PI ) ); // Convert to DEGREES
}

glylib::plane glylib::get_plane_for_ring(int n, GeometryTopology::CoordinateVector r)
{
    int jval[n];
    int l;
    glylib::plane pval; //plane for returning the values of plane
    glylib::vectormag_3D Rj, Rcos, Rsin, R1, R2, R1xR2, avg_coord;
    R1.i = R1.j = R1.k = R2.i = R2.j = R2.k = 0;
    avg_coord.i = avg_coord.j = avg_coord.k = 0;
    for(l = 0; l < n; l++)
    {
        jval[l] = l + 1;//getting j vals for further equations
    }
    for(l = 0; l < n; l++)
    { //for loop for calculating the Rj vals
        Rj.i = r.at(l)->GetX();
        Rj.j = r.at(l)->GetY();
        Rj.k = r.at(l)->GetZ();
        Rsin.i = (Rj.i * sin(2 * gmml::PI_RADIAN * (jval[l] -1) / n));
        Rsin.j = (Rj.j * sin(2 * gmml::PI_RADIAN * (jval[l] -1) / n));
        Rsin.k = (Rj.k * sin(2 * gmml::PI_RADIAN * (jval[l] -1) / n));
        Rcos.i = (Rj.i * cos(2 * gmml::PI_RADIAN * (jval[l] -1) / n));
        Rcos.j = (Rj.j * cos(2 * gmml::PI_RADIAN * (jval[l] -1) / n));
        Rcos.k = (Rj.k * cos(2 * gmml::PI_RADIAN * (jval[l] -1) / n));
        R1.i = R1.i + Rsin.i;
        R1.j = R1.j + Rsin.j;
        R1.k = R1.k + Rsin.k;
        R2.i = R2.i + Rcos.i;
        R2.j = R2.j + Rcos.j;
        R2.k = R2.k + Rcos.k;
    }
    R1xR2 = get_crossprod(R1, R2);//cross prod of Rjs
    //printf("the R1xR2 is %f %f %f %f\n",R1xR2.i,R1xR2.j,R1xR2.k,R1xR2.d);
    //printf("the R1val is %f  %f\n",R1.i,R2.i);
    //printf("the R2val is %f  %f\n",R1.j,R2.j);
    //printf("the R3val is %f  %f\n",R1.k,R2.k);
    pval.A = (R1xR2.i / R1xR2.d);
    pval.B = (R1xR2.j / R1xR2.d);
    pval.C = (R1xR2.k / R1xR2.d);
    for(l = 0; l < n; l++)
    { //for loop for calculating the avg x,y,z coordinates
        avg_coord.i = avg_coord.i + r.at(l)->GetX();
        avg_coord.j = avg_coord.j + r.at(l)->GetY();
        avg_coord.k = avg_coord.k + r.at(l)->GetZ();
    }
    avg_coord.i = (avg_coord.i / n);
    avg_coord.j = (avg_coord.j / n);
    avg_coord.k = (avg_coord.k / n);
    //printf("the avg coords are %f %f %f\n",avg_coord.i,avg_coord.j,avg_coord.k);
    pval.D = -(pval.A*avg_coord.i+pval.B*avg_coord.j+pval.C*avg_coord.k);
    //printf("the dval is %f\n",pval.D);
    return pval;
}

glylib::vectormag_3D glylib::get_crossprod(glylib::vectormag_3D a, glylib::vectormag_3D b)
{
    glylib::vectormag_3D xp;
    xp.i = a.j * b.k - a.k * b.j;
    xp.j = -(a.i * b.k - a.k * b.i);
    xp.k = a.i * b.j - a.j * b.i;
    xp.d = sqrt(xp.i * xp.i + xp.j * xp.j + xp.k * xp.k);
    return xp;
}

double glylib::get_signed_distance_from_point_to_plane(glylib::plane p, GeometryTopology::Coordinate pt)
{
    double sum, sqroot, d;
    //printf("%f %f %f\n",pt.i,pt.j,pt.k);
    //printf("%f %f %f\n",p.A,p.B,p.C);
    //printf("%f %f %f\n",p.A*pt.i,p.B*pt.j,p.C*pt.k);
    sum = (p.A)*(pt.GetX()) + (p.B)*(pt.GetY()) + (p.C)*(pt.GetZ())+p.D;
    //printf(" sum is %f\n",sum);
    sqroot = sqrt(p.A*p.A + p.B*p.B + p.C*p.C);
    d = sum / sqroot;
    //get absolute value
    return d;
}

bool glylib::does_conformation_match(int jloops, int *sortedplanes, int *atoms, int atoms_size, int *list, int required_matches)
{
    int test = count_conformation_matches(jloops, sortedplanes, atoms, atoms_size, list);
    if(test==required_matches)
        return true;
    else
        return false;
}

int glylib::count_conformation_matches(int jloops, int *sortedplanes, int *atoms, int atoms_size, int *list)
{
    int j = 0;
    int match1 = 0;
    while(j < jloops) // maybe sdsize?
    {
        int k = 0;
        int match = 0;
        while(k < atoms_size)
        {
            match=0;
            if(list[sortedplanes[j]*4] == atoms[k]){
                match++;
            }
            if(list[sortedplanes[j]*4+1] == atoms[k+1]){
                match++;
            }
            if(list[sortedplanes[j]*4+2] == atoms[k+2]){
                match++;
            }
            if(list[sortedplanes[j]*4+3] == atoms[k+3]){
                match++;
            }
            if(match==4){
                k = atoms_size;
            }
            k=k+4;
        }
        j++;
        //printf("match is %d\n",match);
        if(match==4)
        {
            match1++;
        }
    }//j<no_planes
    return match1;
}

/************** myfopen(const char,const char) ****************/
FILE* glylib::myfopen(const char *myfilename,const char *myopentype)
{
    FILE *MYFP;
    MYFP=fopen(myfilename,myopentype);
    if(MYFP==NULL){
        printf("Error opening file.\n");
        printf("Expected location is:\n\n");
        printf("\t%s\n\n",myfilename);
        printf("Attempted file open for %s\n",myopentype);
        printf("Exiting.\n");
        exit(1);
    }
    return MYFP;
}

// This doesn't belong here, but Glycan::Monsaccharide is a messy struct.
GeometryTopology::CoordinateVector gmml::GetCycleAtomCoordinates( Glycan::Monosaccharide* mono ) {
    GeometryTopology::CoordinateVector coordinates;
    for( MolecularModeling::AtomVector::iterator it1 = mono->cycle_atoms_.begin(); it1 != mono->cycle_atoms_.end(); ++it1 )
    {
        MolecularModeling::Atom* atom = ( *it1 );
        GeometryTopology::CoordinateVector atom_coordinates = atom->GetCoordinates();
        for( GeometryTopology::CoordinateVector::iterator it2 = atom_coordinates.begin(); it2 != atom_coordinates.end(); ++it2 )
        {
            coordinates.push_back( ( *it2 ) );
        }
    }
    return coordinates;
}
