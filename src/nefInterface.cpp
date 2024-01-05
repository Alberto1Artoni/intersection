#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <iostream>
#include <sstream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;
typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
typedef CGAL::Polyhedron_3<EK> Polyhedron;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::Nef_polyhedron_3<EK> Nef_polyhedron;



std::stringstream createPoly(int nV, int nF, int nc, 
                          double* vx, double* vy, double* vz, int* conn)
{
   // nV: number of verteces
   // nF: number of faces
   // nc: dimension of the connectivity vector
   // vx,vy,vz: coordinate vertex vectors
   // conn: connectivity vector
   double toll  = 0.0000000001;
   double toll2 = 0.00000001;

   std::stringstream  polyOFF;

   polyOFF.precision(8);

   polyOFF << "OFF\n";
   polyOFF << nV << " " << nF << " " << 0 << "\n\n";        // add the two newlines

   // compute diameter in order to have relative dimensions
   double maxX=-1.0, maxY=-1.0, maxZ=-1.0;

// this computes the diameter of the element
   for (int i=0; i < nV; i++)       
   {
       for (int j=(i+1); j < nV; j++)       
       {
           if ( std::abs(vx[i] - vx[j]) > maxX )
           {
               maxX = std::abs(vx[i]-vx[j]);
           } 
           if ( std::abs(vy[i] - vy[j]) > maxY )
           {
               maxY = std::abs(vy[i]-vy[j]);
           } 
           if ( std::abs(vz[i] - vz[j]) > maxZ )
           {
               maxZ = std::abs(vz[i]-vz[j]);
           } 

       }
   }

   double toll2X,toll2Y,toll2Z;

   toll2X = toll2*maxX;
   toll2Y = toll2*maxY;
   toll2Z = toll2*maxZ;

   for (int i=0; i < nV; i++)       
   {
       for (int j=(i+1); j < nV; j++)       // this loop should merge in small squares
       {
           if ((vx[i] - vx[j] < toll2X) && (vx[i] - vx[j] > -toll2X))
           {
               vx[j] = vx[i];
           }
           if ((vy[i] - vy[j] < toll2Y) && (vy[i] - vy[j] > -toll2Y))
           {
               vy[j] = vy[i];
           }
           if ((vz[i] - vz[j] < toll2Z) && (vz[i] - vz[j] > -toll2Z))
           {
               vz[j] = vz[i];
           }
       }
   }

   double tollX,tollY,tollZ;

   tollX = toll*maxX;
   tollY = toll*maxY;
   tollZ = toll*maxZ;

   for (int i=0; i < nV; i++)
   {
  //   polyOFF << vx[i] << " " << vy[i] << " " << vz[i] << "\n";
  //    forse dovrei fare questa cosa qui in fortran?
       if ((vx[i] < tollX) && (vx[i] > -tollX))
       {
           polyOFF << "0" << " ";
       }
       else
       {
           polyOFF << vx[i] << " ";
       }

       if ((vy[i] < tollY) && (vy[i] > -tollY))
       {
           polyOFF << "0" << " ";
       }
       else
       {
           polyOFF << vy[i] << " ";
       }

       if ((vz[i] < tollZ) && (vz[i] > -tollZ))
       {
           polyOFF << "0" << " ";
       }
       else
       {
           polyOFF << vz[i] << " ";
       }

       polyOFF <<  "\n";
   }


// std::cout << "DEBUG:" << std::endl;
// std::cout << conn << std::endl;
// std::cout << conn[0] << std::endl;
// std::cout << conn[conn[0]] << std::endl;

// std::cout << "Stampo conn:" << std::endl;
// for (int i=0; i < nc ; i++)
// {
//     std::cout << conn[i] << " ";
// }

// std::cout << "\n\n\n" << std::endl;


   for (int i=1; i < (nF+1); i++)
   {
        int dimFace = conn[conn[i]];
        polyOFF << dimFace;

        for (int j=conn[i]+1; j < (conn[i]+dimFace+1); j++)
        {
            polyOFF << " " << (conn[j] - 1) ;     // conversion from Fortan to C
        }
        polyOFF << "\n";
   }

// std::cout << "Checking if conversion was successfull:" << std::endl;
// std::cout << polyOFF.str() << std::endl;

   return polyOFF;

   // qua costruisco il nef polyhedron: meglio la stringa
}

// @todo: create a routine that takes the OFF thing are returns the polyhedral structur
// copy the fortran routine
// mmm: per ora non lo faccio, non capisco quanto mi faccia guadagnare
  
void convertPoly( Surface_mesh poly, 
                  int& nV, int& nF, int& nc,
                  double** vx, double** vy, double** vz, int** conn)
{

    nV = poly.number_of_vertices();
    nF = poly.number_of_faces();
    nc = 0;

    // allocate double
 // double* vx = (double*) malloc(nV * sizeof(double));
 // double* vy = (double*) malloc(nV * sizeof(double));
 // double* vz = (double*) malloc(nV * sizeof(double));
 //
    *vx = (double*) malloc(nV * sizeof(double));
    *vy = (double*) malloc(nV * sizeof(double));
    *vz = (double*) malloc(nV * sizeof(double));

    for(vertex_descriptor vd : poly.vertices())
    {
        (*vx)[vd] = CGAL::to_double(poly.point(vd).x());
        (*vy)[vd] = CGAL::to_double(poly.point(vd).y());
        (*vz)[vd] = CGAL::to_double(poly.point(vd).z());
    }

    nc = nF + nF;

    int* faceSizes = (int*) malloc(nF * sizeof(int));

    // get nc dimensions
    for(face_descriptor fd : poly.faces())
    {
        int faceSize = 0;
        Surface_mesh::Halfedge_index hf = poly.halfedge(fd);
        for(Surface_mesh::Halfedge_index hi : halfedges_around_face(hf, poly))      // dunno how to get the face size
        {
            faceSize = faceSize + 1;
        }
        faceSizes[fd] = faceSize;
        nc = nc + faceSize;
    }

    // allocate the nc vector
    // int* conn = (int*) malloc((nc+1) * sizeof(int));
      
       
    *conn = (int*) malloc((nc+1) * sizeof(int));

    (*conn)[0] = nF;
    (*conn)[1] = nF+1;
    for (int i=2; i<=nF; i++)
    {
        (*conn)[i] = faceSizes[i-2] + (*conn)[i-1] + 1;
//      std::cout << i << " " << (*conn)[i] << std::endl;
    }

    for(face_descriptor fd : poly.faces())
    {
        Surface_mesh::Halfedge_index hf = poly.halfedge(fd);
        int counter = 0;
        for(Surface_mesh::Halfedge_index hi : halfedges_around_face(hf, poly))
        {
          Surface_mesh::Vertex_index vi = target(hi, poly);
          counter = counter + 1;
          int vertexLabel(vi);
          (*conn)[(*conn)[fd+1] + counter] = vi;
//        std::cout << (*conn)[fd+1]+counter << " " << (*conn)[(*conn)[fd+1] + counter] << std::endl;
        }
        (*conn)[(*conn)[fd+1]] = counter;   // fd starts from zero
//      std::cout << std::endl;
    }

    free(faceSizes);

//  std::cout << nc <<  " conn: " ;
//  for (int i=0; i<=nc; i++)
//  {
//      std::cout << (*conn)[i] << " ";
//  }
//  std::cout << std::endl;


}
      

extern "C" void   nefInterface(
                          int* nVAcu, int* nFAcu, int* ncAcu, 
                          double* vxAcu, double* vyAcu, double* vzAcu, int* connAcu,
                          int* nVFlu, int* nFFlu, int* ncFlu, 
                          double* vxFlu, double* vyFlu, double* vzFlu, int* connFlu,
                          int* nV, int* nF, int* nc,
                          double** vx, double** vy, double** vz, int** conn)
{


    Polyhedron polyAcu, polyFlu;
    
    std::stringstream ssAcu;
    ssAcu = createPoly(*nVAcu, *nFAcu, *ncAcu, vxAcu, vyAcu, vzAcu, connAcu);
//  std::cout << "building polyhedron Acoustic element inside nef" << std::endl;
    ssAcu >> polyAcu;

//  std::cout << "call ssAcu" << std::endl;
//  std::cout << ssAcu.str() << std::endl;


    std::stringstream ssFlu;
    ssFlu = createPoly(*nVFlu, *nFFlu, *ncFlu, vxFlu, vyFlu, vzFlu, connFlu);
//  std::cout << "building polyhedron fluid element inside nef" << std::endl;
    ssFlu >> polyFlu;

//  std::cout << "call ssFLu" << std::endl;
//  std::cout << ssFlu.str() << std::endl;

    CGAL::Polygon_mesh_processing::triangulate_faces(polyAcu);
    Nef_polyhedron nefAcu(polyAcu);

    //  std::cout << "building Fluid element inside nef" << std::endl;
    CGAL::Polygon_mesh_processing::triangulate_faces(polyFlu);
    Nef_polyhedron nefFlu(polyFlu);

    Nef_polyhedron polyIntersection=nefAcu * nefFlu;

    Surface_mesh output;
    CGAL::convert_nef_polyhedron_to_polygon_mesh(polyIntersection, output);

//  CGAL::IO::write_OFF(std::cout, output, CGAL::parameters::stream_precision(17));

    if (output.number_of_vertices() > 0)
    {

        convertPoly( output, *nV, *nF, *nc , &(*vx), &(*vy), &(*vz), &(*conn));

    }
    else
    {
//      std::cout << "intersection not happening" << std::endl;
        *nV = -1;
    }


}
