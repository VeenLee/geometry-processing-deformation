//#include "biharmonic_precompute.h"
//#include "biharmonic_solve.h"
//#include "arap_precompute.h"
//#include "arap_single_iteration.h"
//#include <igl/min_quad_with_fixed.h>
//#include <igl/read_triangle_mesh.h>
//#include <igl/opengl/glfw/Viewer.h>
//#include <igl/project.h>
//#include <igl/unproject.h>
//#include <igl/snap_points.h>
//#include <igl/unproject_onto_mesh.h>
//#include <Eigen/Core>
//#include <iostream>
//#include <stack>
//
//// Undoable
//struct State
//{
//  // Rest and transformed control points
//  Eigen::MatrixXd CV, CU;
//  bool placing_handles = true;
//} s;
//
//int main(int argc, char *argv[])
//{
//  // Undo Management
//  std::stack<State> undo_stack,redo_stack;
//  const auto push_undo = [&](State & _s=s)
//  {
//    undo_stack.push(_s);
//    // clear
//    redo_stack = std::stack<State>();
//  };
//  const auto undo = [&]()
//  {
//    if(!undo_stack.empty())
//    {
//      redo_stack.push(s);
//      s = undo_stack.top();
//      undo_stack.pop();
//    }
//  };
//  const auto redo = [&]()
//  {
//    if(!redo_stack.empty())
//    {
//      undo_stack.push(s);
//      s = redo_stack.top();
//      redo_stack.pop();
//    }
//  };
//
//  Eigen::MatrixXd V,U;
//  Eigen::MatrixXi F;
//  long sel = -1;
//  Eigen::RowVector3f last_mouse;
//  igl::min_quad_with_fixed_data<double> biharmonic_data, arap_data;
//  Eigen::SparseMatrix<double> arap_K;
//
//  // Load input meshes
//  igl::read_triangle_mesh(
//    (argc>1?argv[1]:"./data/test.ply"),V,F);
//  U = V;
//  igl::opengl::glfw::Viewer viewer;
//  std::cout<<R"(
//[click]  To place new control point
//[drag]   To move control point
//[space]  Toggle whether placing control points or deforming
//M,m      Switch deformation methods
//U,u      Update deformation (i.e., run another iteration of solver)
//R,r      Reset control points 
//⌘ Z      Undo
//⌘ ⇧ Z    Redo
//)";
//  enum Method
//  {
//    BIHARMONIC = 0,
//    ARAP = 1,
//    NUM_METHODS = 2,
//  } method = BIHARMONIC;
//
//  const auto & update = [&]()
//  {
//    // predefined colors
//    const Eigen::RowVector3d orange(1.0,0.7,0.2);
//    const Eigen::RowVector3d yellow(1.0,0.9,0.2);
//    const Eigen::RowVector3d blue(0.2,0.3,0.8);
//    const Eigen::RowVector3d green(0.2,0.6,0.3);
//    if(s.placing_handles)
//    {
//      viewer.data().set_vertices(V);
//      viewer.data().set_colors(blue);
//      viewer.data().set_points(s.CV,orange);
//    }else
//    {
//      // SOLVE FOR DEFORMATION
//      switch(method)
//      {
//        default:
//        case BIHARMONIC:
//        {
//          Eigen::MatrixXd D;
//          biharmonic_solve(biharmonic_data,s.CU-s.CV,D);
//          U = V+D;
//          break;
//        }
//        case ARAP:
//        {
//          arap_single_iteration(arap_data,arap_K,s.CU,U);
//          break;
//        }
//      }
//      viewer.data().set_vertices(U);
//      viewer.data().set_colors(method==BIHARMONIC?orange:yellow);
//      viewer.data().set_points(s.CU,method==BIHARMONIC?blue:green);
//    }
//    viewer.data().compute_normals();
//  };
//  viewer.callback_mouse_down = 
//    [&](igl::opengl::glfw::Viewer&, int, int)->bool
//  {
//    last_mouse = Eigen::RowVector3f(
//      viewer.current_mouse_x,viewer.core().viewport(3)-viewer.current_mouse_y,0);
//    if(s.placing_handles)
//    {
//      // Find closest point on mesh to mouse position
//      int fid;
//      Eigen::Vector3f bary;
//      if(igl::unproject_onto_mesh(
//        last_mouse.head(2),
//        viewer.core().view,
//        viewer.core().proj, 
//        viewer.core().viewport, 
//        V, F, 
//        fid, bary))
//      {
//        long c;
//        bary.maxCoeff(&c);
//        Eigen::RowVector3d new_c = V.row(F(fid,c));
//        __int64 size = s.CV.size();
//        if(size ==0 || (s.CV.rowwise()-new_c).rowwise().norm().minCoeff() > 0)
//        {
//          push_undo();
//          s.CV.conservativeResize(s.CV.rows()+1,3);
//          // Snap to closest vertex on hit face
//          s.CV.row(s.CV.rows()-1) = new_c;
//          update();
//          return true;
//        }
//      }
//    }else
//    {
//      // Move closest control point
//      Eigen::MatrixXf CP;
//      igl::project(
//        Eigen::MatrixXf(s.CU.cast<float>()),
//        viewer.core().view,
//        viewer.core().proj, viewer.core().viewport, CP);
//      Eigen::VectorXf D = (CP.rowwise()-last_mouse).rowwise().norm();
//      sel = (D.minCoeff(&sel) < 30)?sel:-1;
//      if(sel != -1)
//      {
//        last_mouse(2) = CP(sel,2);
//        push_undo();
//        update();
//        return true;
//      }
//    }
//    return false;
//  };
//
//  viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer &, int,int)->bool
//  {
//    if(sel!=-1)
//    {
//      Eigen::RowVector3f drag_mouse(
//        viewer.current_mouse_x,
//        viewer.core().viewport(3) - viewer.current_mouse_y,
//        last_mouse(2));
//      Eigen::RowVector3f drag_scene,last_scene;
//      igl::unproject(
//        drag_mouse,
//        viewer.core().view,
//        viewer.core().proj,
//        viewer.core().viewport,
//        drag_scene);
//      igl::unproject(
//        last_mouse,
//        viewer.core().view,
//        viewer.core().proj,
//        viewer.core().viewport,
//        last_scene);
//      s.CU.row(sel) += (drag_scene-last_scene).cast<double>();
//      last_mouse = drag_mouse;
//      update();
//      return true;
//    }
//    return false;
//  };
//  viewer.callback_mouse_up = [&](igl::opengl::glfw::Viewer&, int, int)->bool
//  {
//    sel = -1;
//    return false;
//  };
//  viewer.callback_key_pressed = 
//    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
//  {
//    switch(key)
//    {
//      case 'M':
//      case 'm':
//      {
//        method = (Method)(((int)(method)+1)%((int)(NUM_METHODS)));
//        break;
//      }
//      case 'R':
//      case 'r':
//      {
//        push_undo();
//        s.CU = s.CV;
//        break;
//      }
//      case 'U':
//      case 'u':
//      {
//        // Just trigger an update
//        break;
//      }
//      case ' ':
//        push_undo();
//        s.placing_handles ^= 1;
//        if(!s.placing_handles && s.CV.rows()>0)
//        {
//          // Switching to deformation mode
//          s.CU = s.CV;
//          Eigen::VectorXi b;
//          igl::snap_points(s.CV,V,b);
//          // PRECOMPUTATION FOR DEFORMATION
//          biharmonic_precompute(V,F,b,biharmonic_data);
//          arap_precompute(V,F,b,arap_data,arap_K);
//        }
//        break;
//      default:
//        return false;
//    }
//    update();
//    return true;
//  };
//
//  // Special callback for handling undo
//  viewer.callback_key_down = 
//    [&](igl::opengl::glfw::Viewer &, unsigned char key, int mod)->bool
//  {
//    if(key == 'Z' && (mod & GLFW_MOD_SUPER))
//    {
//      (mod & GLFW_MOD_SHIFT) ? redo() : undo();
//      update();
//      return true;
//    }
//    return false;
//  };
//  viewer.callback_pre_draw = 
//    [&](igl::opengl::glfw::Viewer &)->bool
//  {
//    if(viewer.core().is_animating && !s.placing_handles && method == ARAP)
//    {
//      arap_single_iteration(arap_data,arap_K,s.CU,U);
//      update();
//    }
//    return false;
//  };
//  viewer.data().set_mesh(V,F);
//  viewer.data().show_lines = false;
//  viewer.core().is_animating = true;
//  viewer.data().face_based = true;
//  update();
//  viewer.launch();
//  return EXIT_SUCCESS;
//}




//#include <igl/colon.h>
//#include <igl/harmonic.h>
//#include <igl/readOBJ.h>
//#include <igl/readDMAT.h>
//#include <igl/opengl/glfw/Viewer.h>
//#include <algorithm>
//#include <iostream>
//
//double bc_frac = 1.0;
//double bc_dir = -0.03;
////bool deformation_field = false;
//Eigen::MatrixXd V, U, V_bc, U_bc;
//Eigen::VectorXd Z;
//Eigen::MatrixXi F;
//Eigen::VectorXi b;
//
//bool pre_draw(igl::opengl::glfw::Viewer& viewer)
//{
//    using namespace Eigen;
//    // Determine boundary conditions
//    if (viewer.core().is_animating)
//    {
//        bc_frac += bc_dir;
//        bc_dir *= (bc_frac >= 1.0 || bc_frac <= 0.0 ? -1.0 : 1.0);
//    }
//
//    const MatrixXd U_bc_anim = V_bc + bc_frac * (U_bc - V_bc);
//    //if (deformation_field)
//    //{
//        MatrixXd D;
//        MatrixXd D_bc = U_bc_anim - V_bc;
//        igl::harmonic(V, F, b, D_bc, 2, D);
//        U = V + D;
//    //}
//    //else
//    //{
//    //    igl::harmonic(V, F, b, U_bc_anim, 2., U);
//    //}
//    viewer.data().set_vertices(U);
//    viewer.data().compute_normals();
//    return false;
//}
//
//bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mods)
//{
//    switch (key)
//    {
//    case ' ':
//        viewer.core().is_animating = !viewer.core().is_animating;
//        return true;
//    //case 'D':
//    //case 'd':
//    //    deformation_field = !deformation_field;
//    //    return true;
//    }
//    return false;
//}
//
//int main(int argc, char* argv[])
//{
//    using namespace Eigen;
//    using namespace std;
//    //igl::readOBJ("./data/decimated-max.obj", V, F);
//    igl::readPLY("./data/test2.ply", V, F);
//    U = V;
//    // S(i) = j: j<0 (vertex i not in handle), j >= 0 (vertex i in handle j)
//    VectorXi S;
//    //igl::readDMAT("./data/decimated-max-selection.dmat", S);
//    igl::readDMAT("./data/test2.dmat", S);
//    igl::colon<int>(0, V.rows() - 1, b);
//    b.conservativeResize(stable_partition(b.data(), b.data() + b.size(),
//        [&S](int i)->bool {return S(i) >= 0; }) - b.data());
//
//    // Boundary conditions directly on deformed positions
//    U_bc.resize(b.size(), V.cols());
//    V_bc.resize(b.size(), V.cols());
//    for (int bi = 0; bi < b.size(); bi++)
//    {
//        V_bc.row(bi) = V.row(b(bi));
//        switch (S(b(bi)))
//        {
//        case 0:
//            // Don't move handle 0
//            U_bc.row(bi) = V.row(b(bi));
//            break;
//        case 1:
//            // move handle 1 down
//            U_bc.row(bi) = V.row(b(bi)) + RowVector3d(0, -10, 0);
//            break;
//        case 2:
//        default:
//            // move other handles forward
//            U_bc.row(bi) = V.row(b(bi)) + RowVector3d(0, 0, -25);
//            break;
//        }
//    }
//
//    // Pseudo-color based on selection
//    MatrixXd C(F.rows(), 3);
//    RowVector3d purple(80.0 / 255.0, 64.0 / 255.0, 255.0 / 255.0);
//    RowVector3d gold(255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0);
//    for (int f = 0; f < F.rows(); f++)
//    {
//        if (S(F(f, 0)) >= 0 && S(F(f, 1)) >= 0 && S(F(f, 2)) >= 0)
//        {
//            C.row(f) = purple;
//        }
//        else
//        {
//            C.row(f) = gold;
//        }
//    }
//
//    // Plot the mesh with pseudocolors
//    igl::opengl::glfw::Viewer viewer;
//    viewer.data().set_mesh(U, F);
//    viewer.data().show_lines = false;
//    viewer.data().set_colors(C);
//    viewer.core().trackball_angle = Eigen::Quaternionf(sqrt(2.0), 0, sqrt(2.0), 0);
//    viewer.core().trackball_angle.normalize();
//    viewer.callback_pre_draw = &pre_draw;
//    viewer.callback_key_down = &key_down;
//    //viewer.core().is_animating = true;
//    viewer.core().animation_max_fps = 30.;
//    cout <<
//        "Press [space] to toggle deformation." << endl /*<<
//        "Press 'd' to toggle between biharmonic surface or displacements." << endl*/;
//    viewer.launch();
//}






#include <igl/colon.h>
#include <igl/directed_edge_orientations.h>
#include <igl/directed_edge_parents.h>
#include <igl/forward_kinematics.h>
#include <igl/PI.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/dqs.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/arap.h>
#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>


typedef
std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> >
RotationList;

const Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
Eigen::MatrixXd V, U;
Eigen::MatrixXi F;
Eigen::VectorXi S, b;
Eigen::RowVector3d mid;
double anim_t = 0.0;
double anim_t_dir = 0.03;
igl::ARAPData arap_data;

bool pre_draw(igl::opengl::glfw::Viewer& viewer)
{
    using namespace Eigen;
    using namespace std;
    MatrixXd bc(b.size(), V.cols());
    for (int i = 0; i < b.size(); i++)
    {
        bc.row(i) = V.row(b(i));
        switch (S(b(i)))
        {
        case 0:
        {
            //head
            const double r = mid(0) * 0.25;
            bc(i, 0) += r * sin(0.5 * anim_t * 2. * igl::PI);
            bc(i, 1) -= r + r * cos(igl::PI + 0.5 * anim_t * 2. * igl::PI);
            break;
        }
        case 1:
        {
            //foot 1
            const double r = mid(1) * 0.15;
            bc(i, 1) += r + r * cos(igl::PI + 0.15 * anim_t * 2. * igl::PI);
            bc(i, 2) -= r * sin(0.15 * anim_t * 2. * igl::PI);
            break;
        }
        //case 2:
        //{
        //    //foot 2
        //    const double r = mid(1) * 0.15;
        //    bc(i, 2) += r + r * cos(igl::PI + 0.35 * anim_t * 2. * igl::PI);
        //    bc(i, 0) += r * sin(0.35 * anim_t * 2. * igl::PI);
        //    break;
        //}
        default:
            break;
        }
    }
    igl::arap_solve(bc, arap_data, U);
    viewer.data().set_vertices(U);
    viewer.data().compute_normals();
    if (viewer.core().is_animating)
    {
        anim_t += anim_t_dir;
    }
    return false;
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mods)
{
    switch (key)
    {
    case ' ':
        viewer.core().is_animating = !viewer.core().is_animating;
        return true;
    }
    return false;
}

int main(int argc, char* argv[])
{
    using namespace Eigen;
    using namespace std;
    //igl::readOFF("./data/decimated-knight.off", V, F);
    //igl::readPLY("./data/test2.ply", V, F);
    igl::readPLY("./data/test.ply", V, F);
    U = V;
    //igl::readDMAT("./data/decimated-knight-selection.dmat", S);
    //igl::readDMAT("./data/test2.dmat", S);
    igl::readDMAT("./data/test.dmat", S);

    // vertices in selection
    //该函数在封闭区间[low,high]内生成size等距的值，step为步长
    igl::colon<int>(0, V.rows() - 1, b);
    //分为两组，满足predicate的排列在前，Return value:Iterator to the first element of the second group
    auto S1 = stable_partition(b.data(), b.data() + b.size(), [](int i)->bool {return S(i) >= 0; });
    auto bSize = S1 - b.data();
    b.conservativeResize(bSize);
    // Centroid
    mid = 0.5 * (V.colwise().maxCoeff() + V.colwise().minCoeff());
    // Precomputation
    arap_data.max_iter = 100;
    igl::arap_precomputation(V, F, V.cols(), b, arap_data);

    // Set color based on selection
    MatrixXd C(F.rows(), 3);
    RowVector3d purple(80.0 / 255.0, 64.0 / 255.0, 255.0 / 255.0);
    RowVector3d gold(255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0);
    for (int f = 0; f < F.rows(); f++)
    {
        if (S(F(f, 0)) >= 0 && S(F(f, 1)) >= 0 && S(F(f, 2)) >= 0)
        {
            C.row(f) = purple;
        }
        else
        {
            C.row(f) = gold;
        }
    }

    // Plot the mesh with pseudocolors
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(U, F);
    viewer.data().set_colors(C);
    viewer.callback_pre_draw = &pre_draw;
    viewer.callback_key_down = &key_down;
    viewer.core().is_animating = false;
    viewer.core().animation_max_fps = 30.;
    cout <<
        "Press [space] to toggle animation" << endl;
    viewer.launch();
}