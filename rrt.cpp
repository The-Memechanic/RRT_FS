#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
#include <SFML/Graphics.hpp>
#include "geometry.h"

#include "tracks/track20.cpp"
#include "tracks/track1.cpp"
#include "tracks/FSG23.cpp"

using namespace std;

const int WIDTH = 256;
const int HEIGHT = 144;
const int FINISH_WIDTH = 8;             // Width of the 'finish line' rectangle
const int FINISH_HEIGHT = 5;            // Height of the 'finish line' rectangle

////////////////////////////////////////
//      Best overall parameters       //
//        JUMP_SIZE = [2...4]         //
//       DISK_SIZE = * [2...5]        //
//   RECENT_NODES_WINDOW = [3...10]   //
////////////////////////////////////////

////////////////////////////////////////
//      Blazing fast parameters       //
//            (DANGEROUS)             //
//          JUMP_SIZE = 5+            //
//          DISK_SIZE = *1            //
//      RECENT_NODES_WINDOW = 1       //
////////////////////////////////////////

double JUMP_SIZE = 2;                 // Maximum distance to jump towards a random point (larger values lead to faster exploration but less optimized paths)
int DISK_SIZE_MULTIPLIER = 3;         // Multiplier for the disk size, which is the radius around a node to search for nearby nodes to rewire
double DISK_SIZE = JUMP_SIZE * DISK_SIZE_MULTIPLIER;     // Circle radius around which we fetch nearby nodes to rewire (larger values lead to more optimized paths, but more execution time and may start going backwards)

const double CAR_WIDTH = 1.2;               // Width of the car for collision checks (resulting path is the midway path if car width is almost equal to the track width)
const double CAR_LENGTH = 1.35;             // Length of the car for collision checks (if the length is too small, the car may face the wall and get stuck)

const int EXTRA_ITERATIONS = 500;           // Extra iterations to keep exploring post success, possibly leading to better paths

int RECENT_NODES_WINDOW = 5;          // Number of furthest nodes to consider expanding, leads to more exploration but also more execution time
                                            // If this value is too low, the algorithm may get stuck facing a wall

int iterations = 0;                         // Number of iterations

Point start;

vector<Point> nodes;                // Nodes of the RRT tree
vector<int> currentPathNodes;       // Nodes of the current best path are stored here
vector<int> parent, nearby;         // Parent of each node, and nearby nodes for rewiring
vector<double> cost, jumps;         // Cost of each node, and jump distance to random point
int nodeCnt = 0, goalIndex = -1;

sf::CircleShape startingPoint;
bool pathFound = false;

Track track;
Point startDirection = {1.0, 0.0};  // Facing right initially I assume
vector<Point> nodeDirection;        // Stores direction of each node

vector<Point> bluePoints, yellowPoints, car_start_points;

void initializeTrack() {
    // Initialize the track with blue and yellow cones
    for (const auto& p : bluePoints) track.addLeftPoint(p);
    for (const auto& p : yellowPoints) track.addRightPoint(p);

    // Start and stop inside the track
    start = car_start_points[0];
}

void prepareInput() {
    // Make starting and ending point circles ready
    startingPoint.setRadius(1);
    startingPoint.setFillColor(sf::Color(208, 0, 240));
    startingPoint.setPosition(start.x - startingPoint.getRadius(), start.y - startingPoint.getRadius());
}

void draw(sf::RenderWindow& window) {
    sf::Vertex line[2];
    sf::CircleShape nodeCircle;

    // Draw track boundaries
    sf::VertexArray leftBoundary(sf::LineStrip, track.leftBoundary.size() + 1);
    for (int i = 0; i < track.leftBoundary.size(); i++) {
        leftBoundary[i].position = sf::Vector2f(track.leftBoundary[i].x, track.leftBoundary[i].y);
        leftBoundary[i].color = sf::Color::Green;
    }
    leftBoundary[track.leftBoundary.size()] = leftBoundary[0]; // close loop
    window.draw(leftBoundary);

    sf::VertexArray rightBoundary(sf::LineStrip, track.rightBoundary.size() + 1);
    for(int i = 0; i < track.rightBoundary.size(); i++) {
        rightBoundary[i].position = sf::Vector2f(track.rightBoundary[i].x, track.rightBoundary[i].y);
        rightBoundary[i].color = sf::Color::Green;
    }
    rightBoundary[track.rightBoundary.size()] = rightBoundary[0]; // close loop
    window.draw(rightBoundary);

    // Draw edges between nodes
    for(int i = (int)nodes.size() - 1; i; i--) {
        Point par = nodes[parent[i]];
        line[0] = sf::Vertex(sf::Vector2f(par.x, par.y));
        line[1] = sf::Vertex(sf::Vector2f(nodes[i].x, nodes[i].y));
        window.draw(line, 2, sf::Lines);
    }

    window.draw(startingPoint);

    // If destination is reached then path is retraced and drawn
    if(pathFound) {
        int node = goalIndex;
        while(parent[node] != node) {
            int par = parent[node];
            line[0] = sf::Vertex(sf::Vector2f(nodes[par].x, nodes[par].y));
            line[1] = sf::Vertex(sf::Vector2f(nodes[node].x, nodes[node].y));
            line[0].color = line[1].color = sf::Color::Red;
            window.draw(line, 2, sf::Lines);
            node = par;
        }
    }

    sf::RectangleShape checkRect;
    sf::Vector2f position = startingPoint.getPosition();
    float width = FINISH_WIDTH * 2;
    float height = FINISH_HEIGHT * 2;
    checkRect.setSize(sf::Vector2f(width, height));
    checkRect.setOrigin(width/2, height/2);  // Center the rectangle
    checkRect.setPosition(position.x, position.y);
    
    // Make it semi-transparent for visibility
    checkRect.setFillColor(sf::Color(255, 255, 0, 50));  // Yellow with transparency
    checkRect.setOutlineColor(sf::Color::Yellow);
    
    window.draw(checkRect);
}

template <typename T> // Returns a random number in [low, high]
T randomCoordinate(T low, T high){
    static random_device random_device;
    static mt19937 engine{random_device()};
    uniform_real_distribution<double> dist(low, high);
    return dist(engine);
}

bool isEdgeObstacleFree(Point a, Point b) {
    Point dir = (b - a).normalized();
    Point side(-dir.y, dir.x); // Perpendicular to direction

    // Car footprint corners centered at a and b
    Point frontA = a + dir * (CAR_LENGTH / 2.0);
    Point rearA  = a - dir * (CAR_LENGTH / 2.0);
    Point frontB = b + dir * (CAR_LENGTH / 2.0);
    Point rearB  = b - dir * (CAR_LENGTH / 2.0);

    // Rectangle corners for swept area (8-point convex hull is overkill; just use quad)
    vector<Point> rect = {
        rearA + side * (CAR_WIDTH / 2.0),
        rearA - side * (CAR_WIDTH / 2.0),
        frontB - side * (CAR_WIDTH / 2.0),
        frontB + side * (CAR_WIDTH / 2.0)
    };

    // Check if any edge of this rectangle intersects track boundaries
    auto intersectsTrack = [&](const vector<Point>& boundary) {
        for (int i = 1; i < boundary.size(); ++i) {
            Point segA = boundary[i - 1];
            Point segB = boundary[i];
            for (int j = 0; j < 4; ++j) {
                if (check_intersection(rect[j], rect[(j + 1) % 4], segA, segB))
                    return true;
            }
        }
        // Close loop
        Point segA = boundary.back();
        Point segB = boundary.front();
        for (int j = 0; j < 4; ++j) {
            if (check_intersection(rect[j], rect[(j + 1) % 4], segA, segB))
                return true;
        }
        return false;
    };

    return !intersectsTrack(track.leftBoundary) && !intersectsTrack(track.rightBoundary);
}


void extractCurrentPath() {
    currentPathNodes.clear();
    int idx = goalIndex;
    //cout << "Extracted Path: ";
    while (true) {
        currentPathNodes.push_back(idx);
        //cout << "(" << nodes[idx].x << ", " << nodes[idx].y << ")";
        if (idx == parent[idx]) break;
        //cout << " <- ";
        idx = parent[idx];
    }
    //cout << endl;
}

void printPathAsVector() {
    if (!pathFound || goalIndex == -1) return;

    vector<Point> pathPoints;
    int current = goalIndex;
    
    // Build the vector with the path points
    while (true) {
        pathPoints.push_back(nodes[current]);
        if (current == parent[current]) break;  // Base case
        current = parent[current];
    }

    std::reverse(pathPoints.begin(), pathPoints.end());

    std::cout << "PATH_VECTOR: vector<Point> pathPoints = {";
    std::cout << "{" << pathPoints[0].x << ", " << pathPoints[0].y << "}";
    for (size_t i = 1; i < pathPoints.size(); ++i) {
        std::cout << ", {" << pathPoints[i].x << ", " << pathPoints[i].y << "}";
    }
    std::cout << "};";
}

// Math...
static double angleBetween(const Point& v1, const Point& v2) {
    double dot = v1.x * v2.x + v1.y * v2.y;
    double det = v1.x * v2.y - v1.y * v2.x;
    double angle = std::atan2(det, dot) * 180.0 / M_PI;
    return std::abs(angle);
}

double calculateMaxAngleChange() {
    if (!pathFound || goalIndex == -1) return 0.0;

    vector<Point> pathPoints;
    int current = goalIndex;
    
    // Build the vector with the path points
    while (true) {
        pathPoints.push_back(nodes[current]);
        if (current == parent[current]) break;  // Base case
        current = parent[current];
    }

    // Check this just to make sure
    if (pathPoints.size() < 3) return 0.0;

    double maxAngleChange = 0.0;

    // Calculate angles between consecutive segments
    for (size_t i = 1; i < currentPathNodes.size() - 1; i++) {
        Point prev = pathPoints[i-1];
        Point curr = pathPoints[i];
        Point next = pathPoints[i+1];

        Point vec1 = (curr - prev);  // Previous segment vector
        Point vec2 = (next - curr);  // Current segment vector

        double angle = angleBetween(vec1, vec2);
        maxAngleChange = std::max(maxAngleChange, angle);
    }

    return maxAngleChange;
}

double calculatePathLengthFromNode(int nodeIndex) {
    if (nodeIndex < 0 || nodeIndex >= nodes.size()) return 0.0;

    double totalLength = 0.0;
    int current = nodeIndex;
    int previous = -1;  // Track previous node to calculate segment distance

    // Traverse from the given node to the root
    while (true) {
        if (previous != -1) {
            // Add distance from current to previous node
            totalLength += distance(nodes[current], nodes[previous]);
        }
        
        // Stop if we've reached the root node
        if (current == parent[current]) {
            break;
        }
        
        // Move to next node in path
        previous = current;
        current = parent[current];
    }

    return totalLength;
}

double calculatePathLength() {
    return calculatePathLengthFromNode(goalIndex);
}

void checkDestinationReached() {
    sf::Vector2f position = startingPoint.getPosition();
    if (nodes.back().x <= position.x + FINISH_WIDTH && nodes.back().x >= position.x - FINISH_WIDTH &&
        nodes.back().y <= position.y + FINISH_HEIGHT && nodes.back().y >= position.y - FINISH_HEIGHT) { // hardcoded rectangle check seems to work best
        if (cost[nodes.size() - 1] < 20)
            return; // Ignore if the path is too short
        pathFound = true;
        goalIndex = nodeCnt - 1;
        // cout << "At " << iterations << " iterations:" << endl;
        // cout << "Found first path with a distance of " << calculatePathLength() << " units. " << endl << endl;
    }
}

double directionScore(Point from, Point to, Point direction) {
    Point movementVec = to - from;
    double dot = movementVec.dot(direction);
    double magProduct = movementVec.magnitude() * direction.magnitude();
    if (magProduct < EPS) return -1.0; // Bad movement

    return dot / magProduct;  // Cosine: 1 = perfect, 0 = perpendicular, -1 = backwards
}

void simplifyPath() {
    if (!pathFound || goalIndex == -1) return;

    int current = goalIndex; // start from the last node

    while (current != parent[current]) { // parent of path's root is itself
        int bestAncestor = parent[current];
        double bestCost = cost[bestAncestor] + distance(nodes[bestAncestor], nodes[current]);

        // Try to skip over intermediate nodes
        int ancestor = parent[bestAncestor];
        while (ancestor != parent[ancestor]) { // again, this will stop at the root
            if (isEdgeObstacleFree(nodes[ancestor], nodes[current])) {
                double newCost = cost[ancestor] + distance(nodes[ancestor], nodes[current]);

                // Skip only if it truly improves cost
                if (newCost + EPS < bestCost) {
                    bestCost = newCost;
                    bestAncestor = ancestor;
                }
            } else {
                break; // Stop checking on the first boundary break
            }
            ancestor = parent[ancestor]; // next ancestor
        }

        // If we found a better parent, manually rewire the path
        if (bestAncestor != parent[current]) {
            parent[current] = bestAncestor;
            cost[current] = bestCost;
            nodeDirection[current] = (nodes[current] - nodes[bestAncestor]).normalized(); // not really necessary in this part of the algorithm, but let's keep it consistent
        }

        current = parent[current]; // next node in the path (which hopefully has been changed)
    }
}


int randomIndex(int min, int max) {
    static std::random_device rd;
    static std::mt19937 rng(rd());
    std::uniform_int_distribution<int> dist(min, max);
    return dist(rng);
}

void RRT() {
    Point nearestPoint, nextPoint;
    bool updated = false;
    int nearestIndex;
    double minCost;
    nearby.clear();
    jumps.resize(nodeCnt);

    while (!updated) {

        int nodeStart = (pathFound) ? 0 : std::max(0, nodeCnt - RECENT_NODES_WINDOW);
        
        vector<int> candidates;
        for (int i = nodeStart; i < nodeCnt; i++) {
            jumps[i] = randomCoordinate(0.3, 1.0) * JUMP_SIZE;
            candidates.push_back(i);
        }

        if (!candidates.empty()) {
            nearestIndex = candidates[randomIndex(0, candidates.size() - 1)];
            nearestPoint = nodes[nearestIndex];
        }

        nextPoint = stepNear(nearestPoint, nodeDirection[nearestIndex],jumps[nearestIndex]); // go forward

        if (!isEdgeObstacleFree(nearestPoint, nextPoint)) {
            continue;
        }

        // Direction biasing
        double dirScore = directionScore(nearestPoint, nextPoint, nodeDirection[nearestIndex]);

        // Randomly decide to skip this point based on direction score (may not be needed anymore)
        double random = randomCoordinate(0.0, 1.0);
        if (!pathFound && random > dirScore) {
            continue;
        }

        // RRT* rewiring path refinement
        for (int i = 0; i < nodeCnt; i++) {
            if ((nodes[i].distance(nextPoint) - DISK_SIZE) <= EPS &&
                isEdgeObstacleFree(nodes[i], nextPoint) &&
                nodes[i].distance(nodes[0]) >= DISK_SIZE) { // Ignore nodes that are too close to the start
                nearby.push_back(i);
            }
        }

        int par = nearestIndex;
        minCost = cost[par] + distance(nodes[par], nextPoint);

        double costBefore = calculatePathLengthFromNode(par); // Calculate the cost of the path before rewiring

        for (auto nodeIndex : nearby) {
            double c = cost[nodeIndex] + distance(nodes[nodeIndex], nextPoint);
            if ((c - minCost) <= EPS) {
                parent.push_back(nodeIndex);
                cost.push_back(c);
                nodes.push_back(nextPoint);
                nodeDirection.push_back((nextPoint - nodes[nodeIndex]).normalized());

                double newPathLength = calculatePathLengthFromNode(nodeIndex); // Calculate the new path length from this node

                if (newPathLength + EPS < costBefore) {
                    minCost = c;
                    par = nodeIndex; // Rewire to this node
                }

                parent.pop_back();
                cost.pop_back();
                nodes.pop_back();
                nodeDirection.pop_back();
            }
        }

        parent.push_back(par);
        cost.push_back(minCost);
        nodes.push_back(nextPoint);
        nodeDirection.push_back((nextPoint - nodes[par]).normalized());

        nodeCnt++;
        updated = true;

        if (!pathFound) {
            checkDestinationReached();
            if (pathFound) {
                simplifyPath(); // Simplify immediately after finding a path
                extractCurrentPath(); // Save current path for random point generation
            }    
        }
        if (pathFound) {
            simplifyPath(); // Minor improvements, but worth it
            extractCurrentPath(); // We keep extracting the path in case it improves
        }
    }
}


int main(int argc, char* argv[]) {
    bool useWindow = true;
    if (argc > 1 && (std::string(argv[1]) == "-n" || std::string(argv[1]) == "--no-window")) {
        useWindow = false;
    }

    if (argc > 2) {
        JUMP_SIZE = std::stod(argv[2]);
        DISK_SIZE_MULTIPLIER = std::stod(argv[3]);
        DISK_SIZE = JUMP_SIZE * DISK_SIZE_MULTIPLIER;
        RECENT_NODES_WINDOW = std::stoi(argv[4]);
        std::string trackName = argv[5];

        // get the track
        if (trackName == "track20") {
            bluePoints = track20_bluePoints;
            yellowPoints = track20_yellowPoints;
            car_start_points = track20_car_start_points;
        } else if (trackName == "track1") {
            bluePoints = track1_bluePoints;
            yellowPoints = track1_yellowPoints;
            car_start_points = track1_car_start_points;
        } else if (trackName == "FSG23") {
            bluePoints = FSG23_bluePoints;
            yellowPoints = FSG23_yellowPoints;
            car_start_points = FSG23_car_start_points;
        } else {
            std::cerr << "Unknown track: " << trackName << "\n";
            exit(1);
        }
    }

    auto startTime = std::chrono::high_resolution_clock::now();
    bool timeoutOccurred = false;
    int postGoalIterations = 0;

    initializeTrack();
    prepareInput();

    sf::RenderWindow window;
    if (useWindow) {
        window.create(sf::VideoMode::getDesktopMode(), "RRT* Path Planning", sf::Style::Fullscreen);
        sf::View view(sf::FloatRect(-WIDTH / 2.0f, -HEIGHT / 2.0f, WIDTH, HEIGHT));
        window.setView(view);
    }

    nodeCnt = 1;
    nodes.push_back(start);
    nodeDirection.push_back(startDirection);
    parent.push_back(0);
    cost.push_back(0);

    //std::cout << "\nWelcome to RRT* Path Planning\n\n";

    while (!useWindow || window.isOpen()) {

        auto currentTime = std::chrono::high_resolution_clock::now();
        auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count();

        // Check if 2 minutes (120,000 ms) have passed
        if (elapsedMs >= 60000) {
            timeoutOccurred = true;
            break;
        }

        if (useWindow) {
            sf::Event event{};
            while (window.pollEvent(event)) {
                if (event.type == sf::Event::Closed) {
                    window.close();
                    return 0;
                }
            }
        }

        RRT();
        iterations++;

        if (pathFound) {
            postGoalIterations++;
            if (postGoalIterations >= EXTRA_ITERATIONS) {
                break;
            }
        }

        // if (iterations % 300 == 0) {
        //     std::cout << "Iterations: " << iterations << std::endl;
        //     if (!pathFound)
        //         std::cout << "Not reached yet :( \n";
        //     else
        //         std::cout << "Shortest distance till now: " << calculatePathLength() << " units.\n";
        //     std::cout << std::endl;
        // }

        if (useWindow) {
            window.clear();
            draw(window);
            window.display();
        }
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

    // std::cout << "\nPath planning completed!\n";
    // std::cout << "Final path distance: " << calculatePathLength() << " units.\n";
    // std::cout << "Total runtime: " << duration << " milliseconds (" << duration / 1000.0 << " seconds).\n";
    // std::cout << "Total iterations: " << iterations << std::endl;
    // std::cout << "Average time per iteration: " << (duration / iterations) << " milliseconds.\n";

    if (!timeoutOccurred) {
        std::cout << "RESULT,"
            << "JumpSize=" << JUMP_SIZE << ","
            << "DiskMultiplier=" << DISK_SIZE_MULTIPLIER << ","
            << "RecentWindow=" << RECENT_NODES_WINDOW << ","
            << "PathLength=" << calculatePathLength() << ","
            << "MaxAngleDeg=" << calculateMaxAngleChange() << ","
            << "RuntimeMs=" << duration << ","
            << "AvgIterationTimeMs=" << (duration / iterations)
            << std::endl;

        printPathAsVector(); // Print the path as a vector of points
    } else {
        std::cout << "RESULT,"
            << "JumpSize=" << JUMP_SIZE << ","
            << "DiskMultiplier=" << DISK_SIZE_MULTIPLIER << ","
            << "RecentWindow=" << RECENT_NODES_WINDOW << ","
            << "PathLength=" << 0.0 << ","
            << "MaxAngleDeg=" << 0.0 << ","
            << "RuntimeMs=" << 0.0 << ","
            << "AvgIterationTimeMs=" << 0.0
            << std::endl;

        std::cout << "PATH_VECTOR: vector<Point> pathPoints = { {0, 0} };" << std::endl;
    }

    if (useWindow) {
        window.close();
    }

    return 0;
}