/* \author Aaron Brown */
// Create simple 3d highway enviroment using PCL
// for exploring self-driving car sensors

#include "sensors/lidar.h"
#include "render/render.h"
#include "processPointClouds.h"
// using templates for processPointClouds so also include .cpp to help linker
#include "processPointClouds.cpp"
#include <pcl/common/transforms.h>
#include <sstream>
#include <string>
#include "quiz/cluster/kdtree.h"
#include <unordered_set>
#include <chrono>

//2019-11-2
std::string num2str(int i)
{
    std::string res;
    std::stringstream ss;
    ss << i;
    ss >> res;
    return res;
}

//my my my -1
std::unordered_set<int> RansacPlane(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, int maxIterations, float distanceTol)
{
        auto startTime = std::chrono::steady_clock::now();

        std::unordered_set<int> inliersResult;
        srand(time(NULL));

        while(maxIterations--)
        {
            //randomly pick two points
            std::unordered_set<int> inliers;
            while(inliers.size() < 3)
                inliers.insert(rand()%(cloud->points.size()));
            float x1, y1, z1, x2, y2, z2, x3, y3, z3;

            auto itr = inliers.begin();
            x1 = cloud->points[*itr].x;
            y1 = cloud->points[*itr].y;
            z1 = cloud->points[*itr].z;
            itr++;
            x2 = cloud->points[*itr].x;
            y2 = cloud->points[*itr].y;
            z2 = cloud->points[*itr].z;
            itr++;
            x3 = cloud->points[*itr].x;
            y3 = cloud->points[*itr].y;
            z3 = cloud->points[*itr].z;

            float a = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1);
            float b = (z2-z1)*(x3-x1)-(x2-x1)*(z3-z1);
            float c = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
            float d = -(a*x1+b*y1+c*z1);

            for(int index = 0; index < cloud->points.size();index++)
            {
                if(inliers.count(index)>0)
                    continue;

                pcl::PointXYZI point = cloud->points[index];
                float x4 = point.x;
                float y4 = point.y;
                float z4 = point.z;

                float dist = fabs(a*x4+b*y4+c*z4+d)/sqrt(a*a+b*b+c*c);

                if(dist <= distanceTol)
                    inliers.insert(index);
            }

            if(inliers.size()>inliersResult.size())
            {
                inliersResult = inliers;
            }
        }

        auto endTime = std::chrono::steady_clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        std::cout << "Ransac took" << elapsedTime.count() << "milliseconds" << std::endl;

        return inliersResult;
}


void clusterHelper(int indice, const std::vector<std::vector<float>>& points, std::vector<int>& cluster, std::vector<bool>& processed, KdTree* tree, float distanceTol)
{
    processed[indice] = true;
    cluster.push_back(indice);

    std::vector<int> nearest = tree->search(points[indice],distanceTol);
    for(int id : nearest)
    {
        if(!processed[id])
            clusterHelper(id, points, cluster, processed, tree, distanceTol);
    }
}

std::vector<std::vector<int>> euclideanCluster(const std::vector<std::vector<float>>& points, KdTree* tree, float distanceTol)
{

        // TODO: Fill out this function to return list of indices for each clus$

        std::vector<std::vector<int>> clusters;

        std::vector<bool> processed(points.size(), false);

        int i = 0;
        while(i < points.size())
        {
            if(processed[i])
            {
                i++;
                continue;
            }
            std::vector<int> cluster;
            clusterHelper(i, points, cluster, processed, tree, distanceTol);
            clusters.push_back(cluster);
            i++;
        }
        return clusters;
}
//my my my 0

std::vector<Car> initHighway(bool renderScene, pcl::visualization::PCLVisualizer::Ptr& viewer)
{

    Car egoCar( Vect3(0,0,0), Vect3(4,2,2), Color(0,1,0), "egoCar");
    Car car1( Vect3(15,0,0), Vect3(4,2,2), Color(0,0,1), "car1");
    Car car2( Vect3(8,-4,0), Vect3(4,2,2), Color(0,0,1), "car2");	
    Car car3( Vect3(-12,4,0), Vect3(4,2,2), Color(0,0,1), "car3");
  
    std::vector<Car> cars;
    cars.push_back(egoCar);
    cars.push_back(car1);
    cars.push_back(car2);
    cars.push_back(car3);

    if(renderScene)
    {
        renderHighway(viewer);
        egoCar.render(viewer);
        car1.render(viewer);
        car2.render(viewer);
        car3.render(viewer);
    }

    return cars;
}


//cityBlock
//void cityBlock(pcl::visualization::PCLVisualizer::Ptr& viewer)
void cityBlock(pcl::visualization::PCLVisualizer::Ptr& viewer, ProcessPointClouds<pcl::PointXYZI> pointProcessor, pcl::PointCloud<pcl::PointXYZI>::Ptr inputCloud)
{
    //ProcessPointClouds<pcl::PointXYZI> pointProcessor;
    //pcl::PointCloud<pcl::PointXYZI>::Ptr inputCloud = pointProcessor.loadPcd("../src/sensors/data/pcd/data_1/0000000000.pcd");
    inputCloud = pointProcessor.FilterCloud(inputCloud, 0.3, Eigen::Vector4f(-10, -5, -2, 1), Eigen::Vector4f(30, 8, 1, 1));
    //////std::pair<pcl::PointCloud<pcl::PointXYZI>::Ptr, pcl::PointCloud<pcl::PointXYZI>::Ptr> segmentCloud = pointProcessor.SegmentPlane(inputCloud, 25, 0.3);
    //my my my1
    std::unordered_set<int> inliers = RansacPlane(inputCloud, 100, .3);

    pcl::PointCloud<pcl::PointXYZI>::Ptr cloudInliers(new pcl::PointCloud<pcl::PointXYZI>());
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloudOutliers(new pcl::PointCloud<pcl::PointXYZI>());

    for(int index = 0; index < inputCloud->points.size(); index++)
    {
        pcl::PointXYZI point = inputCloud->points[index];
        if(inliers.count(index))
            cloudInliers->points.push_back(point);
        else
            cloudOutliers->points.push_back(point);
    }
    //my my my2

    //renderPointCloud(viewer, segmentCloud.first, "obstCloud", Color(1,0,0));
    //////renderPointCloud(viewer, segmentCloud.second, "planeCloud", Color(0,1,0));
    //my my my 3
    renderPointCloud(viewer, cloudInliers, "planeCloud", Color(0,1,0));
    //my my my 4

    //////std::vector<pcl::PointCloud<pcl::PointXYZI>::Ptr> cloudClusters = pointProcessor.Clustering(segmentCloud.first, 0.53, 10, 500);
    //my my my 5
    std::vector<std::vector<float>> points;
    for(int i=0; i < cloudOutliers->points.size(); i++)
        points.push_back({cloudOutliers->points[i].x, cloudOutliers->points[i].y});
    KdTree* tree = new KdTree;
    for (int i=0; i<points.size(); i++) 
        tree->insert(points[i],i);
    std::vector<std::vector<int>> clusters = euclideanCluster(points, tree, 0.53);
    std::vector<pcl::PointCloud<pcl::PointXYZI>::Ptr> cloudClusters;
    for(std::vector<int> cluster : clusters)
    {
        pcl::PointCloud<pcl::PointXYZI>::Ptr cloudCluster(new pcl::PointCloud<pcl::PointXYZI>());
        for(int indice : cluster)
            cloudCluster->points.push_back(cloudOutliers->points[indice]);
        cloudCluster->width = cloudCluster->points.size();
        cloudCluster->height = 1;
        cloudCluster->is_dense = true;

        cloudClusters.push_back(cloudCluster);
    }
    //my my my 6

    int clusterId = 0;
    std::vector<Color> colors = {Color(1,0,0), Color(1,1,0), Color(0,0,1)};
    for(pcl::PointCloud<pcl::PointXYZI>::Ptr cluster : cloudClusters)
    {
        std::cout << "cluster size";
        pointProcessor.numPoints(cluster);
        renderPointCloud(viewer, cluster, "obstCloud"+std::to_string(clusterId), colors[clusterId%colors.size()]);

        Box box = pointProcessor.BoundingBox(cluster);
        renderBox(viewer, box, clusterId);
        ++clusterId;
    }

    //renderPointCloud(viewer, inputCloud, "cloud");
}


void simpleHighway(pcl::visualization::PCLVisualizer::Ptr& viewer)
{
    // ----------------------------------------------------
    // -----Open 3D viewer and display simple highway -----
    // ----------------------------------------------------
    
    // RENDER OPTIONS
    bool renderScene = false;
    std::vector<Car> cars = initHighway(renderScene, viewer);
    
    // TODO:: Create lidar sensor 
    Lidar* lidar = new Lidar(cars, 0);
    pcl::PointCloud<pcl::PointXYZ>::Ptr inputCloud = lidar->scan();
    //renderRays(viewer, lidar->position, inputCloud);
    renderPointCloud(viewer, inputCloud, "inputCloud");
    // TODO:: Create point processor
    ProcessPointClouds<pcl::PointXYZ> pointProcessor;
    std::pair<pcl::PointCloud<pcl::PointXYZ>::Ptr, pcl::PointCloud<pcl::PointXYZ>::Ptr> segmentCloud = pointProcessor.SegmentPlane(inputCloud, 100, 0.2);
    renderPointCloud(viewer,segmentCloud.first, "obstCloud", Color(1,0,0));
    //renderPointCloud(viewer,segmentCloud.second, "planeCloud", Color(0,1,0));

    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> cloudClusters = pointProcessor.Clustering(segmentCloud.first, 1.0, 3, 30);

    int clusterId = 0;
    std::vector<Color> colors = {Color(1,0,0), Color(1,1,0), Color(0,0,1)};
    for(pcl::PointCloud<pcl::PointXYZ>::Ptr cluster : cloudClusters)
    {
        std::cout << "cluster size " << std::endl;
        pointProcessor.numPoints(cluster);
        //renderPointCloud(viewer, cluster, "obstCloud"+std::to_string(clusterId), colors[clusterId%colors.size()]);

        Box box = pointProcessor.BoundingBox(cluster);//2019-11-2
        renderBox(viewer, box, clusterId);//2019-11-2

        // renderPca(viewer, cluster);//2019-11-2,unuseful
        
        /*
        //  11111111111111111111111111111
        //compute principal directions
        Eigen::Vector4f pcaCentroid;
        pcl::compute3DCentroid(*cluster, pcaCentroid);
        Eigen::Matrix3f covariance;
        computeCovarianceMatrixNormalized(*cluster, pcaCentroid, covariance);
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigen_solver(covariance, Eigen::ComputeEigenvectors);
        Eigen::Matrix3f eigenVectorsPCA = eigen_solver.eigenvectors();
        eigenVectorsPCA.col(2) = eigenVectorsPCA.col(0).cross(eigenVectorsPCA.col(1));

        // Transform the original cloud to the origin
        Eigen::Matrix4f projectionTransform(Eigen::Matrix4f::Identity());
        projectionTransform.block<3,3>(0,0) = eigenVectorsPCA.transpose();
        projectionTransform.block<3,1>(0,3) = -1.f * (projectionTransform.block<3,3>(0,0) * pcaCentroid.head<3>());
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPointsProjected(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::transformPointCloud(*cluster, *cloudPointsProjected, projectionTransform);
        // Get the min and max points of the transformed cloud
        pcl::PointXYZ minPoint, maxPoint;
        pcl::getMinMax3D(*cloudPointsProjected, minPoint, maxPoint);
        const Eigen::Vector3f meanDiagonal = 0.5f*(maxPoint.getVector3fMap() + minPoint.getVector3fMap());

        // final transform
        const Eigen::Quaternionf bboxQuaternion(eigenVectorsPCA);
        const Eigen::Vector3f bboxTransform = eigenVectorsPCA * meanDiagonal + pcaCentroid.head<3>();

        //view
        //int mesh_vp_1, mesh_vp_2, mesh_vp_3, mesh_vp_4;
        //viewer->createViewPort(0.0, 0.5, 0.5, 1.0, mesh_vp_1);
        //viewer->createViewPort(0.5, 0.5, 1.0, 1.0, mesh_vp_2);
        //viewer->createViewPort(0.0, 0, 0.5, 0.5, mesh_vp_3);
        //viewer->createViewPort(0.5, 0, 1.0, 0.5, mesh_vp_4);
        // viewer->addPointCloud(cluster, ColorHandlerXYZ(cluster, 30, 144, 255), "bboxedCloud", mesh_vp_3);
        viewer->addCube(bboxTransform, bboxQuaternion, maxPoint.x - minPoint.x, maxPoint.y - minPoint.y, maxPoint.z - minPoint.z, "bbox"+num2str(clusterId));//mesh_vp_3

        //  22222222222222222222222222222
        */
        
        ++clusterId;
        
    }
    //renderPointCloud(viewer, segmentCloud.second, "planeCloud");
    //renderPointCloud(viewer,segmentCloud.first, "obstCloud", Color(1,0,0));

}


//setAngle: SWITCH CAMERA ANGLE {XY, TopDown, Side, FPS}
void initCamera(CameraAngle setAngle, pcl::visualization::PCLVisualizer::Ptr& viewer)
{

    viewer->setBackgroundColor (0, 0, 0);
    
    // set camera position and angle
    viewer->initCameraParameters();
    // distance away in meters
    int distance = 16;
    
    switch(setAngle)
    {
        case XY : viewer->setCameraPosition(-distance, -distance, distance, 1, 1, 0); break;
        case TopDown : viewer->setCameraPosition(0, 0, distance, 1, 0, 1); break;
        case Side : viewer->setCameraPosition(0, -distance, 0, 0, 0, 1); break;
        case FPS : viewer->setCameraPosition(-10, 0, 0, 0, 0, 1);
    }

    if(setAngle!=FPS)
        viewer->addCoordinateSystem (1.0);
}


int main (int argc, char** argv)
{
    std::cout << "starting enviroment" << std::endl;

    pcl::visualization::PCLVisualizer::Ptr viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    //pcl::visualization::PCLVisualizer::Ptr viewer2 (new pcl::visualization::PCLVisualizer (argc, argv, "PlyViewer"));
    CameraAngle setAngle = FPS;
    initCamera(setAngle, viewer);//2019-11-2
    //simpleHighway(viewer);//2019-11-2
    //cityBlock(viewer);

    //create point cloud processor
    ProcessPointClouds<pcl::PointXYZI> pointProcessorI;

    std::vector<boost::filesystem::path> stream = pointProcessorI.streamPcd("../src/sensors/data/pcd/data_1");
    auto streamIterator = stream.begin();

    pcl::PointCloud<pcl::PointXYZI>::Ptr inputCloudI;

    while (!viewer->wasStopped ())
    {
        //clear viewer
        viewer->removeAllPointClouds();
        viewer->removeAllShapes();

        // load pcd and run obstacle detection process
        inputCloudI = pointProcessorI.loadPcd((*streamIterator).string());
        cityBlock(viewer, pointProcessorI, inputCloudI);

        streamIterator++;
        if(streamIterator == stream.end())
            streamIterator = stream.begin();

        viewer->spinOnce ();
    }
}
