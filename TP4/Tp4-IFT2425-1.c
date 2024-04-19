//------------------------------------------------------
// module  : Tp4-IFT2425-1.c
// author  :
// date    :
// version : 1.0
// language: C++
// note    :
//------------------------------------------------------
//

//------------------------------------------------
// FICHIERS INCLUS -------------------------------
//------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <new>
#include <unistd.h>

/************************************************************************/
/* WINDOWS						          	*/
/************************************************************************/
#include <X11/Xutil.h>

Display   *display;
int	  screen_num;
int 	  depth;
Window	  root;
Visual*	  visual;
GC	  gc;

//------------------------------------------------
// DEFINITIONS -----------------------------------
//------------------------------------------------
#define CARRE(X) ((X)*(X))

#define OUTPUT_FILE "Tp4-Img-I.pgm"
#define VIEW_PGM    "xv"
#define DEBUG 0

//-Cst-Modele
#define X_1 0.0
#define Y_1 1.0
#define X_2 -1.0/sqrt(2.0)
#define Y_2 -1.0/2.0
#define X_3 +1.0/2*sqrt(2.0)
#define Y_3 -1.0/2.0
#define C 0.25
#define R 0.1
#define D 0.3

#define X_1_INI 0.2
#define X_2_INI 0.0
#define X_3_INI -1.6
#define X_4_INI 0.0

//-Cst-Runge-Kutta
#define H            0.0001
#define T_0          0.0
#define T_F         30.0
#define NB_INTERV (T_F-T_0)/H

//-Cst-Image
#define WIDTH  512
#define HEIGHT 512
#define MAX_X  4.0
#define MAX_Y  4.0
#define EVOL_GRAPH 3000

#define WHITE     255
#define GREYWHITE 230
#define GREY      200
#define GREYDARK  120
#define BLACK       0

//------------------------------------------------
// GLOBAL CST ------------------------------------
//------------------------------------------------
float Xmin=0.0;
float Xmax=0.0;
float Ymin=0.0;
float Ymax=0.0;

float xx_1=((WIDTH/MAX_X)*X_1)+(WIDTH/2);
float yy_1=(-(HEIGHT/MAX_Y)*Y_1)+(HEIGHT/2);
float xx_2=((WIDTH/MAX_X)*X_2)+(WIDTH/2);
float yy_2=(-(HEIGHT/MAX_Y)*Y_2)+(HEIGHT/2);
float xx_3=((WIDTH/MAX_X)*X_3)+(WIDTH/2);
float yy_3=(-(HEIGHT/MAX_Y)*Y_3)+(HEIGHT/2);

/************************************************************************/
/* OPEN_DISPLAY()							*/
/************************************************************************/
int open_display()
{
    if ((display=XOpenDisplay(NULL))==NULL)
    { printf("Connection impossible\n");
        return(-1); }

    else
    { screen_num=DefaultScreen(display);
        visual=DefaultVisual(display,screen_num);
        depth=DefaultDepth(display,screen_num);
        root=RootWindow(display,screen_num);
        return 0; }
}

/************************************************************************/
/* FABRIQUE_WINDOW()							*/
/* Cette fonction crée une fenetre X et l'affiche à l'écran.	        */
/************************************************************************/
Window fabrique_window(char *nom_fen,int x,int y,int width,int height,int zoom)
{
    Window                 win;
    XSizeHints      size_hints;
    XWMHints          wm_hints;
    XClassHint     class_hints;
    XTextProperty  windowName, iconName;

    char *name=nom_fen;

    if(zoom<0) { width/=-zoom; height/=-zoom; }
    if(zoom>0) { width*=zoom;  height*=zoom;  }

    win=XCreateSimpleWindow(display,root,x,y,width,height,1,0,255);

    size_hints.flags=PPosition|PSize|PMinSize;
    size_hints.min_width=width;
    size_hints.min_height=height;

    XStringListToTextProperty(&name,1,&windowName);
    XStringListToTextProperty(&name,1,&iconName);
    wm_hints.initial_state=NormalState;
    wm_hints.input=True;
    wm_hints.flags=StateHint|InputHint;
    class_hints.res_name=nom_fen;
    class_hints.res_class=nom_fen;

    XSetWMProperties(display,win,&windowName,&iconName,
                     NULL,0,&size_hints,&wm_hints,&class_hints);

    gc=XCreateGC(display,win,0,NULL);

    XSelectInput(display,win,ExposureMask|KeyPressMask|ButtonPressMask|
                             ButtonReleaseMask|ButtonMotionMask|PointerMotionHintMask|
                             StructureNotifyMask);

    XMapWindow(display,win);
    return(win);
}

/****************************************************************************/
/* CREE_XIMAGE()							    */
/* Crée une XImage à partir d'un tableau de float                          */
/* L'image peut subir un zoom.						    */
/****************************************************************************/
XImage* cree_Ximage(float** mat,int z,int length,int width)
{
    int lgth,wdth,lig,col,zoom_col,zoom_lig;
    float somme;
    unsigned char	 pix;
    unsigned char* dat;
    XImage* imageX;

    /*Zoom positiv*/
    /*------------*/
    if (z>0)
    {
        lgth=length*z;
        wdth=width*z;

        dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
        if (dat==NULL)
        { printf("Impossible d'allouer de la memoire.");
            exit(-1); }

        for(lig=0;lig<lgth;lig=lig+z) for(col=0;col<wdth;col=col+z)
            {
                pix=(unsigned char)mat[lig/z][col/z];
                for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
                    {
                        dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+0)]=pix;
                        dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+1)]=pix;
                        dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+2)]=pix;
                        dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+3)]=pix;
                    }
            }
    } /*--------------------------------------------------------*/

        /*Zoom negatifv*/
        /*------------*/
    else
    {
        z=-z;
        lgth=(length/z);
        wdth=(width/z);

        dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
        if (dat==NULL)
        { printf("Impossible d'allouer de la memoire.");
            exit(-1); }

        for(lig=0;lig<(lgth*z);lig=lig+z) for(col=0;col<(wdth*z);col=col+z)
            {
                somme=0.0;
                for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
                        somme+=mat[lig+zoom_lig][col+zoom_col];

                somme/=(z*z);
                dat[((lig/z)*wdth*4)+((4*(col/z))+0)]=(unsigned char)somme;
                dat[((lig/z)*wdth*4)+((4*(col/z))+1)]=(unsigned char)somme;
                dat[((lig/z)*wdth*4)+((4*(col/z))+2)]=(unsigned char)somme;
                dat[((lig/z)*wdth*4)+((4*(col/z))+3)]=(unsigned char)somme;
            }
    } /*--------------------------------------------------------*/

    imageX=XCreateImage(display,visual,depth,ZPixmap,0,(char*)dat,wdth,lgth,16,wdth*4);
    return (imageX);
}

//------------------------------------------------
// FUNCTIONS -------------------------------------
//------------------------------------------------
//-------------------------//
//-- Matrice de Double ----//
//-------------------------//
//---------------------------------------------------------
// Alloue de la memoire pour une matrice 1d de float
//----------------------------------------------------------
float* dmatrix_allocate_1d(int hsize)
{
    float* matrix;
    matrix=new float[hsize]; return matrix; }

//----------------------------------------------------------
// Alloue de la memoire pour une matrice 2d de float
//----------------------------------------------------------
float** dmatrix_allocate_2d(int vsize,int hsize)
{
    float** matrix;
    float *imptr;

    matrix=new float*[vsize];
    imptr=new float[(hsize)*(vsize)];
    for(int i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
    return matrix;
}

//----------------------------------------------------------
// Libere la memoire de la matrice 1d de float
//----------------------------------------------------------
void free_dmatrix_1d(float* pmat)
{ delete[] pmat; }

//----------------------------------------------------------
// Libere la memoire de la matrice 2d de float
//----------------------------------------------------------
void free_dmatrix_2d(float** pmat)
{ delete[] (pmat[0]);
    delete[] pmat;}

//----------------------------------------------------------
// SaveImagePgm
//----------------------------------------------------------
void SaveImagePgm(char* name,float** mat,int lgth,int wdth)
{
    int i,j;
    char buff[300];
    FILE* fic;

    //--extension--
    strcpy(buff,name);

    //--ouverture fichier--
    fic=fopen(buff,"wb");
    if (fic==NULL)
    { printf("Probleme dans la sauvegarde de %s",buff);
        exit(-1); }
    printf("\n Sauvegarde de %s au format pgm\n",buff);

    //--sauvegarde de l'entete--
    fprintf(fic,"P5");
    fprintf(fic,"\n# IMG Module");
    fprintf(fic,"\n%d %d",wdth,lgth);
    fprintf(fic,"\n255\n");

    //--enregistrement--
    for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
            fprintf(fic,"%c",(char)mat[i][j]);

    //--fermeture fichier--
    fclose(fic);
}

//------------------------------------------------------------------------
// plot_point
//
// Affiche entre x dans [-MAX_X/2  MAX_X/2]
//               y dans [-MAX_Y/2  MAX_Y/2]
//------------------------------------------------------------------------
void plot_point(float** MatPts,float** MatPict,int NbPts)
{
    int x_co,y_co;
    int i,j,k;

    //Init
    for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++)  MatPict[i][j]=GREYWHITE;

    for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++)
        { if ((fabs(i-yy_1)+fabs(j-xx_1))<10) MatPict[i][j]=GREYDARK;
            if ((fabs(i-yy_2)+fabs(j-xx_2))<10) MatPict[i][j]=GREYDARK;
            if ((fabs(i-yy_3)+fabs(j-xx_3))<10) MatPict[i][j]=GREYDARK; }

    //Loop
    for(k=0;k<NbPts;k++)
    { x_co=(int)((WIDTH/MAX_X)*MatPts[k][0]);
        y_co=-(int)((HEIGHT/MAX_Y)*MatPts[k][1]);
        y_co+=(HEIGHT/2);
        x_co+=(WIDTH/2);
        if (DEBUG) printf("[%d::%d]",x_co,y_co);
        if ((x_co<WIDTH)&&(y_co<HEIGHT)&&(x_co>0)&&(y_co>0))
            MatPict[y_co][x_co]=BLACK;
    }
}

//------------------------------------------------------------------------
// Fill_Pict
//------------------------------------------------------------------------
void Fill_Pict(float** MatPts,float** MatPict,int PtsNumber,int NbPts)
{
    int i,j;
    int x_co,y_co;
    int k,k_Init,k_End;

    //Init
    for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++)
        { if (MatPict[i][j]!=GREYWHITE) MatPict[i][j]=GREY;
            if ((fabs(i-yy_1)+fabs(j-xx_1))<10) MatPict[i][j]=GREYDARK;
            if ((fabs(i-yy_2)+fabs(j-xx_2))<10) MatPict[i][j]=GREYDARK;
            if ((fabs(i-yy_3)+fabs(j-xx_3))<10) MatPict[i][j]=GREYDARK; }

    //Loop
    k_Init=PtsNumber;
    k_End=(k_Init+EVOL_GRAPH)%NbPts;
    for(k=k_Init;k<k_End;k++)
    { k=(k%NbPts);
        x_co=(int)((WIDTH/MAX_X)*MatPts[k][0]);
        y_co=-(int)((HEIGHT/MAX_Y)*MatPts[k][1]);
        y_co+=(HEIGHT/2);
        x_co+=(WIDTH/2);
        if ((x_co<WIDTH)&&(y_co<HEIGHT)&&(x_co>0)&&(y_co>0))
            MatPict[y_co][x_co]=BLACK; }
}


//------------------------------------------------
// FONCTIONS TPs----------------------------------
//------------------------------------------------
double magnetic_force_x(double pos_x, double mag_x, double pos_y, double mag_y) {
    return (mag_x - pos_x) / pow(sqrt(pow(mag_x - pos_x, 2) + pow(mag_y - pos_y, 2) + pow(D, 2)), 3);
}

double magnetic_force_y(double pos_x, double mag_x, double pos_y, double mag_y) {
    return (mag_y - pos_y) / pow(sqrt(pow(mag_x - pos_x, 2) + pow(mag_y - pos_y, 2) + pow(D, 2)), 3);
}

void derivatives(double t, double y[], double dydt[]) {
    double fx1 = magnetic_force_x(y[0], X_1, y[1], Y_1);
    double fx2 = magnetic_force_x(y[0], X_2, y[1], Y_2);
    double fx3 = magnetic_force_x(y[0], X_3, y[1], Y_3);

    double fy1 = magnetic_force_x(y[1], Y_1, y[0], X_1);
    double fy2 = magnetic_force_x(y[1], Y_2, y[0], X_2);
    double fy3 = magnetic_force_x(y[1], Y_3, y[0], X_3);

    dydt[0] = y[2]; // dx/dt = vx
    dydt[1] = y[3]; // dy/dt = vy
    dydt[2] = -R * y[2] + fx1 + fx2 + fx3 - C; // dvx/dt
    dydt[3] = -R * y[3] + fy1 + fy2 + fy3 - C; // dvy/dt
}

void rkf45(double y[], double t0, double tf, int num_points, float** MatPts) {
    double h = H; // Start with the maximum step size
    double t = t0;
    double dydt[4], k1[4], k2[4], k3[4], k4[4], k5[4], k6[4], y_temp[4];

    int i, j;
    for (i = 0; i < num_points && t < tf; i++) {
        derivatives(t, y, dydt);
        for (j = 0; j < 4; j++) {
            k1[j] = h * dydt[j];
            y_temp[j] = y[j] + 0.25 * k1[j];
        }

        derivatives(t + 0.25 * h, y_temp, dydt);
        for (j = 0; j < 4; j++) {
            k2[j] = h * dydt[j];
            y_temp[j] = y[j] + (3.0/32.0 * k1[j] + 9.0/32.0 * k2[j]);
        }

        derivatives(t + 3.0/8.0 * h, y_temp, dydt);
        for (j = 0; j < 4; j++) {
            k3[j] = h * dydt[j];
            y_temp[j] = y[j] + (1932.0/2197.0 * k1[j] - 7200.0/2197.0 * k2[j] + 7296.0/2197.0 * k3[j]);
        }

        derivatives(t + 12.0/13.0 * h, y_temp, dydt);
        for (j = 0; j < 4; j++) {
            k4[j] = h * dydt[j];
            y_temp[j] = y[j] + (439.0/216.0 * k1[j] - 8.0 * k2[j] + 3680.0/513.0 * k3[j] - 845.0/4104.0 * k4[j]);
        }

        derivatives(t + h, y_temp, dydt);
        for (j = 0; j < 4; j++) {
            k5[j] = h * dydt[j];
            y_temp[j] = y[j] - (8.0/27.0 * k1[j] + 2.0 * k2[j] - 3544.0/2565.0 * k3[j] + 1859.0/4104.0 * k4[j] - 11.0/40.0 * k5[j]);
        }

        derivatives(t + 0.5 * h, y_temp, dydt);
        for (j = 0; j < 4; j++) {
            k6[j] = h * dydt[j];
            y_temp[j] = y[j] + (16.0/135.0 * k1[j] + 6656.0/12825.0 * k3[j] + 28561.0/56430.0 * k4[j] - 9.0/50.0 * k5[j] + 2.0/55.0 * k6[j]);
        }

        // Update the state only with the 5th order estimation
        for (j = 0; j < 4; j++) {
            y[j] = y_temp[j];
        }

        // Store results
        MatPts[i][0] = y[0]; // x coordinate
        MatPts[i][1] = y[1]; // y coordinate

        t += h;  // Advance time by step size h
    }
}


//----------------------------------------------------------
//----------------------------------------------------------
// MAIN
//----------------------------------------------------------
//----------------------------------------------------------
int main (int argc, char **argv)
{
    int i,j,k;
    int flag_graph;
    int zoom;

    XEvent ev;
    Window win_ppicture;
    XImage *x_ppicture;
    char   nomfen_ppicture[100];
    char BufSystVisu[100];

    //>AllocMemory
    float** MatPict=dmatrix_allocate_2d(HEIGHT,WIDTH);
    float** MatPts=dmatrix_allocate_2d((int)(NB_INTERV),2);

    //>Init
    for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) MatPict[i][j]=GREYWHITE;
    for(i=0;i<2;i++) for(j=0;j<(int)(NB_INTERV);j++) MatPts[i][j]=0.0;
    flag_graph=1;
    zoom=1;


    //---------------------------------------------------------------------
    //>Question 1
    //---------------------------------------------------------------------

    //Il faut travailler ici ...et dans > // FONCTIONS TPs

    //Un exemple ou la matrice de points est remplie
    //par une courbe donné par l'équation d'en bas... et non pas par
    //la solution de l'équation différentielle
    double y[] = {0.2,-1.6,0.05,1};

    /*for(k=0;k<(int)(NB_INTERV);k++)
    {
        MatPts[k][0]=(k/(float)(NB_INTERV))*cos((k*0.0001)*3.14159);
        //MatPts[k][1]=(k/(float)(NB_INTERV))*sin((k*0.001)*3.14159);
        //>on peut essayer la ligne d'en bas aussi
        MatPts[k][1]=(k/(float)(NB_INTERV))*sin((k*0.0001)*3.14159);
    }*/

    rkf45(y, T_0, T_F, NB_INTERV, MatPts);

    // Print or process trajectory points
    for (int i = 0; i < NB_INTERV; i++) {
        if(i%1000 == 0) {
            printf("Point %d: (%f, %f)\n", i, MatPts[i][0], MatPts[i][1]);
        }
    }


    //--Fin Question 1-----------------------------------------------------


    //>Affichage des Points dans MatPict
    plot_point(MatPts,MatPict,(int)(NB_INTERV));

    //>Save&Visu de MatPict
    SaveImagePgm((char*)OUTPUT_FILE,MatPict,HEIGHT,WIDTH);
    strcpy(BufSystVisu,VIEW_PGM);
    strcat(BufSystVisu," ");
    strcat(BufSystVisu,OUTPUT_FILE);
    strcat(BufSystVisu," &");
    system(BufSystVisu);

    //>Affiche Statistique
    printf("\n\n Stat:  Xmin=[%.2f] Xmax=[%.2f] Ymin=[%.2f] Ymax=[%.2f]\n",Xmin,Xmax,Ymin,Ymax);


    //--------------------------------------------------------------------------------
    //-------------- visu sous XWINDOW de l'évolution de MatPts ----------------------
    //--------------------------------------------------------------------------------
    if (flag_graph)
    {
        //>Uuverture Session Graphique
        if (open_display()<0) printf(" Impossible d'ouvrir une session graphique");
        sprintf(nomfen_ppicture,"Évolution du Graphe");
        win_ppicture=fabrique_window(nomfen_ppicture,10,10,HEIGHT,WIDTH,zoom);
        x_ppicture=cree_Ximage(MatPict,zoom,HEIGHT,WIDTH);

        printf("\n\n Pour quitter,appuyer sur la barre d'espace");
        fflush(stdout);

        //>Boucle Evolution
        for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) MatPict[i][j]=GREYWHITE;
        for(k=0;;)
        {
            k=((k+EVOL_GRAPH)%(int)(NB_INTERV));
            Fill_Pict(MatPts,MatPict,k,(int)(NB_INTERV));
            XDestroyImage(x_ppicture);
            x_ppicture=cree_Ximage(MatPict,zoom,HEIGHT,WIDTH);
            XPutImage(display,win_ppicture,gc,x_ppicture,0,0,0,0,x_ppicture->width,x_ppicture->height);
            usleep(10000);  //si votre machine est lente mettre un nombre plus petit
        }
    }

    //>Retour
    printf("\n Fini... \n\n\n");
    return 0;
}




