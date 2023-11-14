#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#define STB_IMAGE_IMPLEMENTATION
#pragma warning(disable:4996)
#include <unordered_map>
#include "shader.h"
#include "camera.h"
#include "basic_camera.h"
#include "pointLight.h"
#include "cube.h"
#include "sphere.h"
#include "stb_image.h"
#include <iostream>
#include <ctime>

using namespace std;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);
void drawCube(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 model, float r, float g, float b);
void getCurrentTime(int& hours, int& minutes, int& seconds);
glm::mat4 transform(float tx, float ty, float tz, float sx, float sy, float sz);
// ************************DRAWING ROOM****************************************
void f_drawing_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube);
void f_drawing_sofa(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& sofa_top, Cube& sofa_foam, Cube& sofa_pillow);
void f_drawing_sofa1(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& sofa_top, Cube& sofa_foam, Cube& sofa_pillow, float rotateY, float tx, float ty, float tz);
void f_drawing_table(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_table);
void f_drawing_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& wall, Cube& door, Shader& ourShader);
void f_drawing_tv(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_tv, Cube& drawing_sound_box, Cube& drawing_cupboard);
void f_drawing_mat(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_mat);
void fan(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& b, Cube& c, float x = 0.0f, float y = 0.0f, float z = 0.0f);
void f_drawing_window(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_window);

//  *************************************BED ROOM 1******************************************
void f_bed_room1_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube, Shader& ourShader, Sphere& clock_bell, Cube& clock, Cube& x);
void f_bed_room1_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube);
void f_bed_room1_bed(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& bed_sheet, Cube& bed_pillow, Cube& bed_texture, Cube& blanket_texture);
void f_bed_room1_almari(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& draw_almari);
void f_bed_room1_window(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_bed_room1_window);
void f_bed_room1_bedtable(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_bed_room1_bedtable, Shader& ourShader);
//  *************************************DINING ROOM******************************************
void f_dining_room_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube);
void f_dining_room_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube);
void f_dining_table(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& dining_table);
void f_dining_frize(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_dining_frize);
void f_dining_window(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& kitchen_window);
void f_dining_light(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Shader& ourShader);
void chair(unsigned int& cubeVAO, Shader& lightingShaderWithTexture, glm::mat4 alTogether, Cube& chair_texture, float tx, float ty, float tz, float rotateY);

//  *************************************KITCHEN***********************************************
void f_kitchen_surface(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& kitchen_surface, Cube& kitchen_surface_top, Cube& kitchen_cupboard1, Cube& kitchen_cupboard2, Cube& kitchen_back_texture);
void f_kitchen_stove(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& kitchen_surface);

//  *************************************BED ROOM 2********************************************
void f_bed_room2_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube);
void f_bed_room2_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube);
void f_bed_room2_book(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_dining_frize);
void f_bed_room2_light(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Shader& ourShader);
void f_bed_room2_dressing_table(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& dressing_right, Cube& dressing_bottom, Cube& dressing_mirror);
void f_bed_room2_bed(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& bed_sheet, Cube& bed_pillow, Cube& bed_texture, Cube& blanket_texture);

// **************************************ROOF***************************************************
void f_roof(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_roof);

// ***************************************BATHROOM**********************************************
void f_bathroom_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& bathroom_top, Cube& bathroom_door, Cube& bathroom_tiles);
void f_bathroom_light(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Shader& ourShader);
void f_toilet(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& toilet);
void f_bathroom_bucket(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& water);

unsigned int loadTexture(char const* path, GLenum textureWrappingModeS, GLenum textureWrappingModeT, GLenum textureFilteringModeMin, GLenum textureFilteringModeMax);
void shaderActivate(Shader& shader);

// settings
const unsigned int SCR_WIDTH = 1000;
const unsigned int SCR_HEIGHT = 800;

// modelling transform
float rotateAngle_X = 0.0;
float rotateAngle_Y = 0.0;
float rotateAngle_Z = 0.0;
float rotateAxis_X = 0.0;
float rotateAxis_Y = 0.0;
float rotateAxis_Z = 1.0;
float translate_X = 0.0;
float translate_Y = 0.0;
float translate_Z = 0.0;
float scale_X = 1.0;
float scale_Y = 1.0;
float scale_Z = 1.0;

// front door
float f_door = 0.0f;
float b1_door = 0.0f;
// bathroom door
float max_bathroom_door_translate = (3.0 + 0.1 - 1.75) / 2.0;
float bathroom_door_translate = 0.0f;

// fan
float rotateFan = 0;
float rotateClock = 0.0f;
bool sign = 1;



// camera
Camera camera(glm::vec3(0.0f, 1.1f, 5.2f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

float eyeX = 0.0, eyeY = 1.0, eyeZ = 3.0;
float lookAtX = 0.0, lookAtY = 0.0, lookAtZ = 0.0;
glm::vec3 V = glm::vec3(0.0f, 1.0f, 0.0f);
BasicCamera basic_camera(eyeX, eyeY, eyeZ, lookAtX, lookAtY, lookAtZ, V);


// positions of the point lights
glm::vec3 pointLightPositions[] = {
    glm::vec3(0.17f, 0.4f, -1.75f),
    glm::vec3(0.0f,  1.5f,  0.0f),
    glm::vec3(0.0f,  1000.0f,  0.0f),
    glm::vec3(0.0f,  3.0f,  0.0f)
};

glm::vec3 point_light_positions[] = {
    glm::vec3(1.45f, 1.3f, 0.1f),
    glm::vec3(1.45f, 1.3f, -3.1f),
    glm::vec3(1.6f, 1.3f, -3.1f),
    glm::vec3(1.6f, 1.3f, 0.5f),
    glm::vec3(2.5f + 1.9f, 0.8f, -0.9f)
};

PointLight pointlight1(

    pointLightPositions[0].x, pointLightPositions[0].y, pointLightPositions[0].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,        // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    1       // light number
);
PointLight pointlight2(

    pointLightPositions[1].x, pointLightPositions[1].y, pointLightPositions[1].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    2       // light number
);

PointLight pointlight3(

    pointLightPositions[2].x, pointLightPositions[2].y, pointLightPositions[2].z,  // position
    0.1f, 0.1f, 0.1f,     // ambient
    0.1f, 0.1f, 0.1f,      // diffuse
    0.1f, 0.1f, 0.1f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    3       // light number
);
PointLight pointlight4(

    pointLightPositions[3].x, pointLightPositions[3].y, pointLightPositions[3].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    4       // light number
);
// ******************************DRAWING_ROOM_LIGHT***********************************
PointLight drawing_light(
    point_light_positions[0].x, point_light_positions[0].y, point_light_positions[0].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    5       // light number
);
PointLight bed_room1_light(
    point_light_positions[1].x, point_light_positions[1].y, point_light_positions[1].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    6       // light number
);
PointLight dining_light(
    point_light_positions[2].x, point_light_positions[2].y, point_light_positions[2].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    7       // light number
);
PointLight bed_room2_light(
    point_light_positions[3].x, point_light_positions[3].y, point_light_positions[3].z,  // position
    0.3f, 0.3f, 0.3f,     // ambient
    0.3f, 0.3f, 0.3f,      // diffuse
    0.3f, 0.3f, 0.3f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    8       // light number
);
PointLight bathroom_light(
    point_light_positions[4].x, point_light_positions[4].y, point_light_positions[4].z,  // position
    0.2f, 0.2f, 0.2f,     // ambient
    0.2f, 0.2f, 0.2f,      // diffuse
    0.2f, 0.2f, 0.2f,         // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    9       // light number
);



// light settings
bool onOffPointToggle = true;
bool onOffSpotToggle = true;
bool onOffDirectToggle = true;
bool ambientToggle = true;
bool diffuseToggle = true;
bool specularToggle = true;

//glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
//glm::mat4 view = camera.GetViewMatrix();
glm::mat4 projection;
glm::mat4 view;

string diffuseMapPath;
string specularMapPath;

class Curve
{
public:
    vector<float> cntrlPoints;
    vector <float> coordinates;
    vector <float> normals;
    vector <int> indices;
    vector <float> vertices;
    const double pi = 3.14159265389;
    const int nt = 40;
    const int ntheta = 20;
    Curve(vector<float>& tmp)
    {
        this->cntrlPoints = tmp;
        this->fishVAO = hollowBezier(cntrlPoints.data(), ((unsigned int)cntrlPoints.size() / 3) - 1);
        cout << cntrlPoints.size() << endl;
        cout << coordinates.size() << endl;
        cout << normals.size() << endl;
        cout << indices.size() << endl;
        cout << vertices.size() << endl;
    }
    ~Curve()
    {
        glDeleteVertexArrays(1, &fishVAO);
        glDeleteVertexArrays(1, &bezierVAO);
        glDeleteBuffers(1, &bezierVBO);
        glDeleteBuffers(1, &bezierEBO);
    }
    void draw(Shader& lightingShader, glm::mat4 model)
    {
        /// Fish
        lightingShader.use();
        lightingShader.setMat4("model", model);
        lightingShader.setVec3("material.ambient", glm::vec3(1.0f, 0.6f, 0.0f));
        lightingShader.setVec3("material.diffuse", glm::vec3(1.0f, 0.6f, 0.0f));
        lightingShader.setVec3("material.specular", glm::vec3(1.0f, 1.0f, 1.0f));
        lightingShader.setFloat("material.shininess", 32.0f);

        glBindVertexArray(fishVAO);
        glDrawElements(GL_TRIANGLES,                    // primitive type
            (unsigned int)indices.size(),          // # of indices
            GL_UNSIGNED_INT,                 // data type
            (void*)0);                       // offset to indices

        // unbind VAO
        glBindVertexArray(0);
        /// End Fish
    }
private:
    unsigned int fishVAO;
    unsigned int bezierVAO;
    unsigned int bezierVBO;
    unsigned int bezierEBO;


    unsigned int drawControlPoints()
    {
        unsigned int controlPointVAO;
        unsigned int controlPointVBO;

        glGenVertexArrays(1, &controlPointVAO);
        glGenBuffers(1, &controlPointVBO);

        glBindVertexArray(controlPointVAO);

        glBindBuffer(GL_ARRAY_BUFFER, controlPointVBO);
        glBufferData(GL_ARRAY_BUFFER, (unsigned int)cntrlPoints.size() * sizeof(float), cntrlPoints.data(), GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        return controlPointVAO;
    }

    long long nCr(int n, int r)
    {
        if (r > n / 2)
            r = n - r; // because C(n, r) == C(n, n - r)
        long long ans = 1;
        int i;

        for (i = 1; i <= r; i++)
        {
            ans *= n - r + i;
            ans /= i;
        }

        return ans;
    }
    void BezierCurve(double t, float xy[2], GLfloat ctrlpoints[], int L)
    {
        double y = 0;
        double x = 0;
        t = t > 1.0 ? 1.0 : t;
        for (int i = 0; i < L + 1; i++)
        {
            long long ncr = nCr(L, i);
            double oneMinusTpow = pow(1 - t, double(L - i));
            double tPow = pow(t, double(i));
            double coef = oneMinusTpow * tPow * ncr;
            x += coef * ctrlpoints[i * 3];
            y += coef * ctrlpoints[(i * 3) + 1];

        }
        xy[0] = float(x);
        xy[1] = float(y);
    }
    unsigned int hollowBezier(GLfloat ctrlpoints[], int L)
    {
        int i, j;
        float x, y, z, r;                //current coordinates
        float theta;
        float nx, ny, nz, lengthInv;    // vertex normal


        const float dtheta = 2 * pi / ntheta;        //angular step size

        float t = 0;
        float dt = 1.0 / nt;
        float xy[2];

        for (i = 0; i <= nt; ++i)              //step through y
        {
            BezierCurve(t, xy, ctrlpoints, L);
            r = xy[0];
            y = xy[1];
            theta = 0;
            t += dt;
            lengthInv = 1.0 / r;

            for (j = 0; j <= ntheta; ++j)
            {
                double cosa = cos(theta);
                double sina = sin(theta);
                z = r * cosa;
                x = r * sina;

                coordinates.push_back(x);
                coordinates.push_back(y);
                coordinates.push_back(z);

                // normalized vertex normal (nx, ny, nz)
                // center point of the circle (0,y,0)
                nx = (x - 0) * lengthInv;
                ny = (y - y) * lengthInv;
                nz = (z - 0) * lengthInv;

                normals.push_back(nx);
                normals.push_back(ny);
                normals.push_back(nz);

                theta += dtheta;
            }
        }
        // generate index list of triangles
        // k1--k1+1
        // |  / |
        // | /  |
        // k2--k2+1

        int k1, k2;
        for (int i = 0; i < nt; ++i)
        {
            k1 = i * (ntheta + 1);     // beginning of current stack
            k2 = k1 + ntheta + 1;      // beginning of next stack

            for (int j = 0; j < ntheta; ++j, ++k1, ++k2)
            {
                // k1 => k2 => k1+1
                indices.push_back(k1);
                indices.push_back(k2);
                indices.push_back(k1 + 1);

                // k1+1 => k2 => k2+1
                indices.push_back(k1 + 1);
                indices.push_back(k2);
                indices.push_back(k2 + 1);
            }
        }

        size_t count = coordinates.size();
        for (int i = 0; i < count; i += 3)
        {
            //cout << count << ' ' << i + 2 << endl;
            vertices.push_back(coordinates[i]);
            vertices.push_back(coordinates[i + 1]);
            vertices.push_back(coordinates[i + 2]);

            vertices.push_back(normals[i]);
            vertices.push_back(normals[i + 1]);
            vertices.push_back(normals[i + 2]);
        }

        glGenVertexArrays(1, &bezierVAO);
        glBindVertexArray(bezierVAO);

        // create VBO to copy vertex data to VBO
        glGenBuffers(1, &bezierVBO);
        glBindBuffer(GL_ARRAY_BUFFER, bezierVBO);           // for vertex data
        glBufferData(GL_ARRAY_BUFFER,                   // target
            (unsigned int)vertices.size() * sizeof(float), // data size, # of bytes
            vertices.data(),   // ptr to vertex data
            GL_STATIC_DRAW);                   // usage

        // create EBO to copy index data
        glGenBuffers(1, &bezierEBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bezierEBO);   // for index data
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,           // target
            (unsigned int)indices.size() * sizeof(unsigned int),             // data size, # of bytes
            indices.data(),               // ptr to index data
            GL_STATIC_DRAW);                   // usage

        // activate attrib arrays
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);

        // set attrib arrays with stride and offset
        int stride = 24;     // should be 24 bytes
        glVertexAttribPointer(0, 3, GL_FLOAT, false, stride, (void*)0);
        glVertexAttribPointer(1, 3, GL_FLOAT, false, stride, (void*)(sizeof(float) * 3));

        // unbind VAO, VBO and EBO
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

        return bezierVAO;
    }

};

Curve* bucket;

vector<float>Bucket = {
    -1.0000, 1.8450, 5.1000,
        -0.9600, 1.9100, 5.1000,
        -0.9150, 1.8700, 5.1000,
        -0.9000, 1.7350, 5.1000,
        -0.8600, 1.5950, 5.1000,
        -0.8500, 1.4400, 5.1000,
        -0.8200, 1.2450, 5.1000,
        -0.8100, 1.0600, 5.1000,
        -0.7800, 0.8200, 5.1000,
        -0.7500, 0.6300, 5.1000,
        -0.6900, 0.4650, 5.1000,
        -0.5750, 0.4750, 5.1000,
        -0.4350, 0.4600, 5.1000,
        -0.3100, 0.4450, 5.1000,
        -0.2150, 0.4600, 5.1000,
        -0.1250, 0.4600, 5.1000,
        -0.0350, 0.4600, 5.1000,
        -0.0050, 0.4650, 5.1000,
};

// timing
float deltaTime = 0.0f;    // time between current frame and last frame
float lastFrame = 0.0f;

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Flat Interior", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // build and compile our shader zprogram
    // ------------------------------------
    Shader lightingShader("vertexShaderForPhongShading.vs", "fragmentShaderForPhongShading.fs");
    Shader lightingShaderWithTexture("vertexShaderForPhongShadingWithTexture.vs", "fragmentShaderForPhongShadingWithTexture.fs");
    Shader ourShader("vertexShader.vs", "fragmentShader.fs");

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------

    float cube_vertices[] = {
        // positions      // normals
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        1.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        1.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        0.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,

        1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,

        0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,

        0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 1.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f,

        1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f,
        1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f,

        0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
        1.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f
    };
    unsigned int cube_indices[] = {
        0, 3, 2,
        2, 1, 0,

        4, 5, 7,
        7, 6, 4,

        8, 9, 10,
        10, 11, 8,

        12, 13, 14,
        14, 15, 12,

        16, 17, 18,
        18, 19, 16,

        20, 21, 22,
        22, 23, 20
    };

    unsigned int cubeVAO, cubeVBO, cubeEBO;
    glGenVertexArrays(1, &cubeVAO);
    glGenBuffers(1, &cubeVBO);
    glGenBuffers(1, &cubeEBO);

    glBindVertexArray(cubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cube_vertices), cube_vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cube_indices), cube_indices, GL_STATIC_DRAW);


    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // vertex normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)12);
    glEnableVertexAttribArray(1);

    // second, configure the light's VAO (VBO stays the same; the vertices are the same for the light object which is also a 3D cube)
    unsigned int lightCubeVAO;
    glGenVertexArrays(1, &lightCubeVAO);
    glBindVertexArray(lightCubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    // note that we update the lamp's position attribute's stride to reflect the updated buffer data
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);


    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // drawing_floor--------------------------------------------------------
    diffuseMapPath = "images/floor2.jpg";
    specularMapPath = "images/floor2.jpg";
    Cube drawing_floor = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);

    // fan_center-----------------------------------------------------------
    diffuseMapPath = "images/fan_center.PNG";
    specularMapPath = "images/fan_center.PNG";
    Cube c = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // fan_blade--------------------------------------------------------------------
    diffuseMapPath = "images/f_b.png";
    specularMapPath = "images/f_b.png";
    Cube b = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // drawing_wall-----------------------------------------------------------
    diffuseMapPath = "images/wall_color.PNG";
    specularMapPath = "images/wall_color.PNG";
    Cube drawing_wall = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);

    // drawing_window-----------------------------------------------------------
    diffuseMapPath = "images/drawing_window.jpg";
    specularMapPath = "images/drawing_window.jpg";
    Cube drawing_window = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    //   sofa_top ------------------------------------------------------------
    diffuseMapPath = "images/sofa_top.jpg";
    specularMapPath = "images/sofa_top.jpg";
    Cube sofa_top = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);


    // sofa_foam --------------------------------------------------------------
    diffuseMapPath = "images/sofa_foam.jpg";
    specularMapPath = "images/sofa_foam.jpg";
    Cube sofa_foam = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);

    // sofa_pillow-----------------------------------------------------------------
    diffuseMapPath = "images/sofa_pillow.jpg";
    specularMapPath = "images/sofa_pillow.jpg";
    Cube sofa_pillow = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);

    //  drawing_table---------------------------------------------------------------
    diffuseMapPath = "images/table_top.jpg";
    specularMapPath = "images/table_top.jpg";
    Cube drawing_table = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 6.0f, 4.0f);

    //  drawing_tv-------------------------------------------------------------------
    diffuseMapPath = "images/tv.jpg";
    specularMapPath = "images/tv.jpg";
    Cube drawing_tv = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    //  drawing_sound_box-------------------------------------------------------------
    diffuseMapPath = "images/soundbox.jpg";
    specularMapPath = "images/soundbox.jpg";
    Cube drawing_sound_box = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 4.0f, 1.0f);

    //  drawing_door-------------------------------------------------------------
    diffuseMapPath = "images/door.PNG";
    specularMapPath = "images/door.PNG";
    Cube door = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // drawing_cupboard-----------------------------------------------------------------
    diffuseMapPath = "images/soundbox.jpg";
    specularMapPath = "images/soundbox.jpg";
    Cube drawing_cupboard = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 4.0f, 2.0f);

    // drawing_mat --------------------------------------------------------------------
    diffuseMapPath = "images/mat.jpg";
    specularMapPath = "images/mat.jpg";
    Cube drawing_mat = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // bed_room1_almari --------------------------------------------------------------------
    diffuseMapPath = "images/almari.png";
    specularMapPath = "images/almari.png";
    Cube drawing_almari = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 4.0f);

    // bed_room1_window ---------------------------------------------------------------------
    diffuseMapPath = "images/bed_room1_window.jpeg";
    specularMapPath = "images/bed_room1_window.jpeg";
    Cube drawing_bed_room1_window = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // bed_room1_bedtable ---------------------------------------------------------------------
    diffuseMapPath = "images/bed_side_table.jpg";
    specularMapPath = "images/bed_side_table.jpg";
    Cube drawing_bed_room1_bedtable = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 2.0f);

    // bed_room1_bed_sheet ---------------------------------------------------------------------
    diffuseMapPath = "images/bed_sheet.jpg";
    specularMapPath = "images/bed_sheet.jpg";
    Cube bed_sheet = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 3.0f);
    
    // bed_room1_bed_texture ---------------------------------------------------------------------
    diffuseMapPath = "images/bed_texture.PNG";
    specularMapPath = "images/bed_texture.PNG";
    Cube bed_texture = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 3.0f, 1.0f);
    
    // bed_room1_blanket_texture ---------------------------------------------------------------------
    diffuseMapPath = "images/clock.png";
    specularMapPath = "images/clock.png";
    Cube clock = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // bed_room1_blanket_texture ---------------------------------------------------------------------
    diffuseMapPath = "images/blanket_texture.jpg";
    specularMapPath = "images/blanket_texture.jpg";
    Cube blanket_texture = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 2.0f, 2.0f);

    // bed_room1_bed_pillow ---------------------------------------------------------------------
    diffuseMapPath = "images/bed_pillow.jpg";
    specularMapPath = "images/bed_pillow.jpg";
    Cube bed_pillow = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 3.0f);

    // dining_frize ---------------------------------------------------------------------
    diffuseMapPath = "images/frize.PNG";
    specularMapPath = "images/frize.PNG";
    Cube dining_frize = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    
    // dining_frize ---------------------------------------------------------------------
    diffuseMapPath = "images/dining_window.jpg";
    specularMapPath = "images/dining_window.jpg";
    Cube dining_window = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);


    // kitchen_surface ---------------------------------------------------------------------
    diffuseMapPath = "images/kitchen_surface.jpg";
    specularMapPath = "images/kitchen_surface.jpg";
    Cube kitchen_surface = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 8.0f, 3.0f);

    // kitchen_surface_top ---------------------------------------------------------------------
    diffuseMapPath = "images/kitchen_surface_top.jpg";
    specularMapPath = "images/kitchen_surface_top.jpg";
    Cube kitchen_surface_top = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 5.0f, 2.0f);
    
    // kitchen_stove ---------------------------------------------------------------------
    diffuseMapPath = "images/stove.PNG";
    specularMapPath = "images/stove.PNG";
    Cube kitchen_stove = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // kitchen_back ---------------------------------------------------------------------
    diffuseMapPath = "images/kitchen_back_texture.jpg";
    specularMapPath = "images/kitchen_back_texture.jpg";
    Cube kitchen_back_texture = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // kitchen_cupboard ---------------------------------------------------------------------
    diffuseMapPath = "images/kitchen_cupboard.PNG";
    specularMapPath = "images/kitchen_cupboard.PNG";
    /*diffuseMapPath = "images/kitchen_back_texture.jpg";
    specularMapPath = "images/kitchen_back_texture.jpg";*/
    Cube kitchen_cupboard1 = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 5.0f, 1.0f);
    Cube kitchen_cupboard2 = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 2.0f, 1.0f);

    // bed_room2_bookshelf
    diffuseMapPath = "images/book.jpg";
    specularMapPath = "images/book.jpg";
    Cube bookshelf = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // roof
    diffuseMapPath = "images/roof.jpg";
    specularMapPath = "images/roof.jpg";
    Cube roof = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 12.0f, 8.0f);

    // chair
    diffuseMapPath = "images/chair_texture.PNG";
    specularMapPath = "images/chair_texture.PNG";
    Cube chair_texture = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // dining_table_texture
    diffuseMapPath = "images/dining_table_texture.PNG";
    specularMapPath = "images/dining_table_texture.PNG";
    Cube dining_table_texture = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 2.0f, 1.0f);

    // dressing_right
    diffuseMapPath = "images/dressing_right.PNG";
    specularMapPath = "images/dressing_right.PNG";
    Cube dressing_right = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // dressing_bottom
    diffuseMapPath = "images/dressing_bottom.PNG";
    specularMapPath = "images/dressing_bottom.PNG";
    Cube dressing_bottom = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // dressing_mirror
    diffuseMapPath = "images/mirror.PNG";
    specularMapPath = "images/mirror.PNG";
    Cube dressing_mirror = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // bathroom_top
    diffuseMapPath = "images/bathroom_top.PNG";
    specularMapPath = "images/bathroom_top.PNG";
    Cube bathroom_top = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // bathroom_top
    diffuseMapPath = "images/bathroom_door.PNG";
    specularMapPath = "images/bathroom_door.PNG";
    Cube bathroom_door = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    
    // bathroom_tiles
    diffuseMapPath = "images/bathroom_tiles.PNG";
    specularMapPath = "images/bathroom_tiles.PNG";
    Cube bathroom_tiles = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 4.0f, 2.0f);

    // bathroom_toilet
    diffuseMapPath = "images/toilet.PNG";
    specularMapPath = "images/toilet.PNG";
    Cube bathroom_toilet = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);

    // bathroom_toilet
    diffuseMapPath = "images/water.PNG";
    specularMapPath = "images/water.PNG";
    Cube water = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 2.0f, 1.0f);

    // bathroom_toilet
    diffuseMapPath = "images/dd.PNG";
    specularMapPath = "images/dd.PNG";
    Cube dd = Cube(loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR), 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);


    // drawing_point_light
    Sphere clock_bell = Sphere();

    //ourShader.use();
    //lightingShader.use();

    Curve buck(Bucket);
    bucket = &buck;

    pointlight2.turnOff();

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        // --------------------
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // be sure to activate shader when setting uniforms/drawing objects
        lightingShader.use();
        lightingShader.setVec3("viewPos", camera.Position);

        //// point light 1
        //pointlight1.setUpPointLight(lightingShader);
        //// point light 2
        //pointlight2.setUpPointLight(lightingShader);
        //// point light 3
        //pointlight3.setUpPointLight(lightingShader);
        //// point light 4
        //pointlight4.setUpPointLight(lightingShader);


        // *************************************DRAWING ROOM LIGHT********************************
        //drawing_light.setUpPointLight(lightingShader);

        
        // activate shader
        lightingShader.use();

        // pass projection matrix to shader (note that in this case it could change every frame)
        // glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        //glm::mat4 projection = glm::ortho(-2.0f, +2.0f, -1.5f, +1.5f, 0.1f, 100.0f);
        projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        lightingShader.setMat4("projection", projection);

        // camera/view transformation
        // glm::mat4 view = camera.GetViewMatrix();
        //glm::mat4 view = basic_camera.createViewMatrix();
        view = camera.GetViewMatrix();
        lightingShader.setMat4("view", view);

        // Modelling Transformation
        glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
        glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        lightingShader.setMat4("model", model);

        //glBindVertexArray(cubeVAO);
        //glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        //glDrawArrays(GL_TRIANGLES, 0, 36);

        /*bed(cubeVAO, lightingShader, model,0);
        bed(cubeVAO, lightingShader, model, 2);*/
        //drawingRoom(cubeVAO, lightingShader, model, lightingShaderWithTexture);
        
        /*building(cubeVAO, lightingShader, model);
        road(cubeVAO, lightingShader, model);*/

        // ********************************DRAWING ROOM************************************
        f_drawing_floor(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_floor);
        f_drawing_wall(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_wall, door, ourShader);
        f_drawing_sofa(cubeVAO, lightingShader, model, lightingShaderWithTexture, sofa_top, sofa_foam, sofa_pillow);
        f_drawing_sofa1(cubeVAO, lightingShader, model, lightingShaderWithTexture, sofa_top, sofa_foam, sofa_pillow, 90.0, 0, 0, 0);
        f_drawing_table(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_table);
        f_drawing_tv(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_tv, drawing_sound_box, drawing_cupboard);
        f_drawing_mat(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_mat);
        fan(cubeVAO, lightingShader, model, lightingShaderWithTexture, b, c, 0.0f, 0.0f, 0.0f);
        f_drawing_window(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_window);
        
        //  **********************************BED ROOM 1*************************************
        f_bed_room1_wall(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_wall, ourShader, clock_bell, clock, dd);
        f_bed_room1_floor(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_floor);
        f_bed_room1_bed(cubeVAO, lightingShader, model, lightingShaderWithTexture, bed_sheet, bed_pillow, bed_texture, blanket_texture);
        f_bed_room1_almari(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_almari);
        f_bed_room1_window(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_bed_room1_window);
        f_bed_room1_bedtable(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_bed_room1_bedtable, ourShader);
        fan(cubeVAO, lightingShader, model, lightingShaderWithTexture, b, c, 0.0f, 0.0f, -2.8f);
        
        // ***************************************DINING ROOM***************************************
        f_dining_room_wall(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_wall);
        f_dining_room_floor(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_floor);
        f_dining_table(cubeVAO, lightingShader, model, lightingShaderWithTexture, dining_table_texture);
        f_dining_frize(cubeVAO, lightingShader, model, lightingShaderWithTexture, dining_frize);
        f_dining_window(cubeVAO, lightingShader, model, lightingShaderWithTexture, dining_window);
        fan(cubeVAO, lightingShader, model, lightingShaderWithTexture, b, c, 3.0f, 0.0f, -2.5f);
        f_dining_light(cubeVAO, lightingShader, model, lightingShaderWithTexture, ourShader);
        chair(cubeVAO, lightingShaderWithTexture, model, chair_texture, 9.5f, 1.5f, -7.0f, 0.0f);
        chair(cubeVAO, lightingShaderWithTexture, model, chair_texture, 11.0f, 1.5f, -7.0f, 0.0f);
        
        // ****************************************KITCHEN******************************************
        f_kitchen_surface(cubeVAO, lightingShader, model, lightingShaderWithTexture, kitchen_surface, kitchen_surface_top, kitchen_cupboard1, kitchen_cupboard2, kitchen_back_texture);
        f_kitchen_stove(cubeVAO, lightingShader, model, lightingShaderWithTexture, kitchen_stove);
        
        // *****************************************BED ROOM 2**************************************
        f_bed_room2_wall(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_wall);
        f_bed_room2_floor(cubeVAO, lightingShader, model, lightingShaderWithTexture, drawing_floor);
        f_bed_room2_book(cubeVAO, lightingShader, model, lightingShaderWithTexture, bookshelf);
        fan(cubeVAO, lightingShader, model, lightingShaderWithTexture, b, c, 3.5f, 0.0f, 0.5f);
        f_bed_room2_light(cubeVAO, lightingShader, model, lightingShaderWithTexture, ourShader);
        f_bed_room2_dressing_table(cubeVAO, lightingShader, model, lightingShaderWithTexture, dressing_right, dressing_bottom, dressing_mirror);
        f_bed_room2_bed(cubeVAO, lightingShader, model, lightingShaderWithTexture, bed_sheet, bed_pillow, bed_texture, blanket_texture);
        
        // ********************************************ROOF******************************************
        f_roof(cubeVAO, lightingShader, model, lightingShaderWithTexture, roof);

        // ***********************************************BATHROOM***********************************
        f_bathroom_wall(cubeVAO, lightingShader, model, lightingShaderWithTexture, bathroom_top, bathroom_door, bathroom_tiles);
        f_bathroom_light(cubeVAO, lightingShader, model, lightingShaderWithTexture, ourShader);
        f_toilet(cubeVAO, lightingShader, model, lightingShaderWithTexture, bathroom_toilet);
        f_bathroom_bucket(cubeVAO, lightingShader, model, lightingShaderWithTexture, water);

        // also draw the lamp object(s)
        ourShader.use();
        ourShader.setMat4("projection", projection);
        ourShader.setMat4("view", view);

        // we now draw as many light bulbs as we have point lights.
        glBindVertexArray(lightCubeVAO);
        //for (unsigned int i = 0; i < 4; i++)
        //{
        //    model = glm::mat4(1.0f);
        //    model = glm::translate(model, pointLightPositions[i]);
        //    model = glm::scale(model, glm::vec3(0.2f)); // Make it a smaller cube
        //    ourShader.setMat4("model", model);
        //    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
        //    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        //    //glDrawArrays(GL_TRIANGLES, 0, 36);
        //}

        ////texture
        //glm::mat4 modelTexture = glm::mat4(1.0f);
        //glm::mat4 translate = glm::mat4(1.0f);
        //glm::mat4 scale = glm::mat4(0.5f);

        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X+1, translate_Y+1, translate_Z + 1));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X * 0.5, scale_Y * 0.5, scale_Z * 0.5));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
        
        lightingShaderWithTexture.use();
        lightingShaderWithTexture.setVec3("viewPos", camera.Position);
        lightingShaderWithTexture.setMat4("view", view);
        lightingShaderWithTexture.setMat4("projection", projection);

        lightingShaderWithTexture.use();
        // point light 1
        pointlight1.setUpPointLight(lightingShaderWithTexture);
        // point light 2
        pointlight2.setUpPointLight(lightingShaderWithTexture);
        // point light 3
        pointlight3.setUpPointLight(lightingShaderWithTexture);
        // point light 4
        pointlight4.setUpPointLight(lightingShaderWithTexture);

        //  *****************************DRAWING_LIGHT*****************************
        drawing_light.setUpPointLight(lightingShaderWithTexture);
        //  ******************************BED ROOM 1********************************
        bed_room1_light.setUpPointLight(lightingShaderWithTexture);
        //   ******************************DINING ROOM******************************
        dining_light.setUpPointLight(lightingShaderWithTexture);
        //   *******************************BED ROOM2*******************************
        bed_room2_light.setUpPointLight(lightingShaderWithTexture);
        //   ******************************BATHROOM LIGHT*********************************
        bathroom_light.setUpPointLight(lightingShaderWithTexture);

        /*diffuseMapPath = "images/emoji.png";
        specularMapPath = "images/emoji.png";
        diffMap = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
        specMap = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
        Cube cube = Cube(diffMap, specMap, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
        cube.drawCubeWithTexture(lightingShaderWithTexture, model);*/


        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    glDeleteVertexArrays(1, &cubeVAO);
    glDeleteVertexArrays(1, &lightCubeVAO);
    glDeleteBuffers(1, &cubeVBO);
    glDeleteBuffers(1, &cubeEBO);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

void shaderActivate(Shader& shader)
{
    shader.use();
    shader.setVec3("viewPos", camera.Position);
    shader.setMat4("view", view);
    shader.setMat4("projection", projection);
}


void drawCube(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 model = glm::mat4(1.0f), float r = 1.0f, float g = 1.0f, float b = 1.0f)
{
    lightingShader.use();

    lightingShader.setVec3("material.ambient", glm::vec3(r, g, b));
    lightingShader.setVec3("material.diffuse", glm::vec3(r, g, b));
    lightingShader.setVec3("material.specular", glm::vec3(0.5f, 0.5f, 0.5f));
    lightingShader.setFloat("material.shininess", 32.0f);

    lightingShader.setMat4("model", model);

    glBindVertexArray(cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
}

// *******************DRAWING ROOM*****************************
void f_drawing_sofa(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& sofa_top, Cube& sofa_foam, Cube& sofa_pillow)
{
    float baseHeight = 0.2f;
    float width = 1.2f;
    float length = 0.4f;
    
    // base
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0.0, -0.5));
    translate2 = glm::translate(model, glm::vec3(width-0.5- 0.5, 0.0, -1.0));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.546, 0.335, 0.316);
    shaderActivate(lightingShaderWithTexture);
    sofa_top.drawCubeWithTexture(lightingShaderWithTexture, model);
    
    // foam
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.03, 0.3));
    translate2 = glm::translate(model, glm::vec3(width - 0.19 - 0.5, 0.2, -1.0+0.05));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.03, 0.3));
    translate2 = glm::translate(model, glm::vec3(width - 0.195-0.3 - 0.5, 0.2, -1.0 + 0.05));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.03, 0.3));
    translate2 = glm::translate(model, glm::vec3(width - 0.2-0.3-0.3 - 0.5, 0.2, -1.0 + 0.05));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // back
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, 0.1));
    translate2 = glm::translate(model, glm::vec3(width - 0.5 - 0.5, 0.2, -1.15));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    // right
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.12, baseHeight, 0.4));
    translate2 = glm::translate(model, glm::vec3(width+0.04 - 0.5, 0.2, -1.0));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.12, baseHeight, 0.4));
    translate2 = glm::translate(model, glm::vec3(width - 1.04 - 0.5, 0.2, -1.0));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left pillo
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 0.1, 0.01));
    translate2 = glm::translate(model, glm::vec3(width-0.8 - 0.5, 0.23, -1.09));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    sofa_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);

    // right pillo
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 0.1, 0.01));
    translate2 = glm::translate(model, glm::vec3(width - 0.19 - 0.5, 0.23, -1.09));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    sofa_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

}


void f_drawing_sofa1(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& sofa_top, Cube& sofa_foam, Cube& sofa_pillow, float rotateY, float tx, float ty, float tz)
{
    float baseHeight = 0.2f;
    float width = 1.2f;
    float length = 0.4f;

    // base
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0.0, -0.5));
    translate2 = glm::translate(model, glm::vec3(width - 0.5 - 0.5, 0.0, -1.0));
    glm::mat4 rotateYMatrix = glm::rotate(model, glm::radians(rotateY), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.546, 0.335, 0.316);
    shaderActivate(lightingShaderWithTexture);
    sofa_top.drawCubeWithTexture(lightingShaderWithTexture, model);

    // foam
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.03, 0.3));
    translate2 = glm::translate(model, glm::vec3(width - 0.19 - 0.5, 0.2, -1.0 + 0.05));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.03, 0.3));
    translate2 = glm::translate(model, glm::vec3(width - 0.195 - 0.3 - 0.5, 0.2, -1.0 + 0.05));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.03, 0.3));
    translate2 = glm::translate(model, glm::vec3(width - 0.2 - 0.3 - 0.3 - 0.5, 0.2, -1.0 + 0.05));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // back
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, 0.1));
    translate2 = glm::translate(model, glm::vec3(width - 0.5 - 0.5, 0.2, -1.15));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    // right
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.12, baseHeight, 0.4));
    translate2 = glm::translate(model, glm::vec3(width + 0.04 - 0.5, 0.2, -1.0));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.12, baseHeight, 0.4));
    translate2 = glm::translate(model, glm::vec3(width - 1.04 - 0.5, 0.2, -1.0));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.2, 0.335, 0.5);
    sofa_foam.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left pillo
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 0.1, 0.01));
    translate2 = glm::translate(model, glm::vec3(width - 0.8 - 0.5, 0.23, -1.09));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    sofa_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);

    // right pillo
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 0.1, 0.01));
    translate2 = glm::translate(model, glm::vec3(width - 0.19 - 0.5, 0.23, -1.09));
    model = alTogether * rotateYMatrix * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    sofa_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

}


void f_drawing_table(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_table)
{
    float baseHeight = 0.02;
    float width = 0.7;
    float length = 0.4;

    // base
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0.0, -0.5));
    translate2 = glm::translate(model, glm::vec3(width - 0.5, 0.2, -0.2));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.1, 0.1, 0.1);
    shaderActivate(lightingShaderWithTexture);
    drawing_table.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    // leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.2, 0.05));
    translate2 = glm::translate(model, glm::vec3(width-0.175, 0.0, -0.35));
    model = alTogether * translate2 * scale * translate;
    drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);

    // leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.2, 0.05));
    translate2 = glm::translate(model, glm::vec3(width - 0.175, 0.0, -0.025));
    model = alTogether * translate2 * scale * translate;
    drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);

    // leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.2, 0.05));
    translate2 = glm::translate(model, glm::vec3(width - 0.825, 0.0, -0.35));
    model = alTogether * translate2 * scale * translate;
    drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);

    // leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.2, 0.05));
    translate2 = glm::translate(model, glm::vec3(width - 0.825, 0.0, -0.025));
    model = alTogether * translate2 * scale * translate;
    drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
    
}


void f_drawing_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube &cube)
{
    float baseHeight = 0.01;
    float width = 3.0;
    float length = 3.0;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width+1, baseHeight, length+1));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_drawing_mat(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_mat)
{
    float baseHeight = 0.01;
    float width = 0.6;
    float length = 0.4;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(0.0, 1.0, 0.5));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 1.0, 0.0, 0.0);
    shaderActivate(lightingShaderWithTexture);
    drawing_mat.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_drawing_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube, Cube& door, Shader& ourShader)
{
    float baseHeight = 1.5;
    float width = 0.01;
    float length = 3.0;

    // right
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width+0.1, baseHeight, length+0.1));
    translate = glm::translate(model, glm::vec3(1.5, 0.0, -1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // left
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width+0.1, baseHeight, length+0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // back
    width = 3.0;
    length = 0.01;
    //model = glm::mat4(1.0f);
    //scale = glm::scale(model, glm::vec3(width, baseHeight, length + 0.1));
    //translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    //model = alTogether * translate * scale;
    ////drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //

    // front top
    float baseHeight1 = baseHeight / 3;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight1, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 2*baseHeight1, 1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // front bottom
    float baseHeight2 = baseHeight1 * 2;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width-0.5, baseHeight2, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.0, 0.0, 1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    

    // back top
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight1, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 2 * baseHeight1, -1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // back bottom
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width - 0.5, baseHeight2, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // front door
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.5, baseHeight2, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, 1.5));
    glm::mat4 rotateYMatrix = glm::rotate(model, glm::radians(f_door), glm::vec3(0.0f, 1.0f, 0.0f));
    glm::mat4 translate_origin = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    model = alTogether * translate * rotateYMatrix * translate_origin * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    door.drawCubeWithTexture(lightingShaderWithTexture, model);

    // back door
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.5, baseHeight2, length + 0.1));
    translate = glm::translate(model, glm::vec3(width-2.0, 0.0, -1.5));
    rotateYMatrix = glm::rotate(model, glm::radians(b1_door), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * translate_origin * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    door.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(ourShader);
    
    // point light
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.15, 0.15));
    translate = glm::translate(model, glm::vec3(1.45f, 1.3f, 0.1f));
    model = alTogether * translate * scale;
    ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    shaderActivate(lightingShader);

}


void f_drawing_tv(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_tv, Cube& drawing_sound_box, Cube& drawing_cupboard)
{
    float baseHeight = 0.4;
    float width = 0.6;
    float length = 0.02;

    // tv
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(0.0, 0.5, 1.475));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.0, 0.0, 0.0);
    shaderActivate(lightingShaderWithTexture);
    drawing_tv.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // sound box
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, 0.05, length));
    translate = glm::translate(model, glm::vec3(0.0, 0.4, 1.475));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.0, 0.0, 0.0);
    drawing_sound_box.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // cupboard
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, 0.2, length+0.15));
    translate = glm::translate(model, glm::vec3(0.0, 0.0, 1.475-0.15));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.0, 0.0, 0.0);
    drawing_cupboard.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);
}


void fan(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& b, Cube& c, float x, float y, float z)
{
    float bladel = 1.5;
    float bladew = 0.2;
    float bladeh = 0.01;

    // Center
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 scale2 = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.27, 0.3, 0.27));
    scale2 = glm::scale(model, glm::vec3(0.5, 0.5, 0.5));
    translate = glm::translate(model, glm::vec3(-0.67, 0.0, -0.4));
    translate2 = glm::translate(model, glm::vec3(0.0, 1.35, 0.0));
    translate3 = glm::translate(model, glm::vec3(x, y, z));
    model = alTogether * translate3 * translate2 * scale2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    c.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0.01, 0.0, 0.0));
    glm::mat4 rotateM = glm::rotate(model, glm::radians(45.0f + rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate3 * translate2 * scale2 * rotateM * scale * translate;
    shaderActivate(lightingShaderWithTexture);
    b.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0, 0.01, 0.0));
    rotateM = glm::rotate(model, glm::radians(165.0f + rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate3 * translate2 * scale2 * rotateM * scale * translate;
    b.drawCubeWithTexture(lightingShaderWithTexture, model);


    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0.01, 0.01, 0.0));
    rotateM = glm::rotate(model, glm::radians(285.0f + rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate3 * translate2 * scale2 * rotateM * scale * translate;
    b.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_drawing_window(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_window)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1, 0.7, 0.05));
    glm::mat4 rotateY = glm::rotate(model, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    translate = glm::translate(model, glm::vec3(-0.3, 0.5, -1.4));
    model = alTogether * rotateY * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    drawing_window.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


// *************************BED ROOM 1*************************************
void f_bed_room1_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube, Shader& ourShader, Sphere& clock_bell, Cube& clock, Cube& x)
{
    float baseHeight = 1.5;
    float width = 0.01;
    float length = 3.0;

    // right
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 rotateYMatrix = glm::rotate(model, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(1.5, 0.0, -1.5));
    translate2 = glm::translate(model, glm::vec3(0.0, 0.0, -length));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    // left
    //model = glm::mat4(1.0f);
    //scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight, length + 0.1));
    //translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    //model = alTogether * translate2 * rotateYMatrix * translate * scale;
    ////drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // back
    width = 3.0;
    length = 0.01;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // front top
    baseHeight = baseHeight / 3;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 2 * baseHeight, 1.5));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // front bottom
    baseHeight = baseHeight * 2;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width - 0.5, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.0, 0.0, 1.5));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // point light
    shaderActivate(ourShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.15, 0.15));
    translate = glm::translate(model, glm::vec3(1.45f, 1.3f, -3.1f));
    model = alTogether * translate * scale;
    ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    //shaderActivate(lightingShader);

    // clock
    shaderActivate(lightingShaderWithTexture);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.1, 0.3, 0.3));
    translate = glm::translate(model, glm::vec3(-1.35f , 1.1f + 0.13f - 0.03-0.03, -3.15f + 0.03f- 0.02- 0.05));
    model = alTogether * translate * scale;
    /*ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    */
    //shaderActivate(lightingShader);
    clock.drawCubeWithTexture(lightingShaderWithTexture, model);

    shaderActivate(ourShader);
    if (rotateClock >= 40) {
        sign =0;
    }
    else if (rotateClock <= -40) {
        sign = 1;
    }
    if (sign == 1) {
        rotateClock += 0.5;
    }
    else
        rotateClock -= 0.5;

    glm::mat4 translate_origin = glm::mat4(1.0f);
    glm::mat4 translate4 = glm::mat4(1.0f);
    model = glm::mat4(1.0f);
    translate_origin = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    scale = glm::scale(model, glm::vec3(0.02, -0.3, 0.03));
    translate = glm::translate(model, glm::vec3(-1.35f, 0.83f + 0.3f + 0.1, -3.04f));
    glm::mat4 rotateM = glm::rotate(model, glm::radians(0.0f + rotateClock), glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * translate * rotateM * translate_origin * scale;
    //model = alTogether * translate * translate_origin * rotateM * translate_origin * scale;
    /*ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);*/
    x.drawCubeWithTexture(lightingShaderWithTexture, model);


    shaderActivate(lightingShaderWithTexture);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.02, 0.04, 0.04));
    glm::mat4 rotateb = glm::rotate(model, glm::radians(180.0f), glm::vec3(0.0f, 0.0f, 1.0f));
    translate = glm::translate(model, glm::vec3(-1.35f, 0.8f+0.3f + 0.1, -3.03f));
    translate4 = glm::translate(model, glm::vec3(0, 0.3, 0));
    model = alTogether * translate * rotateM * rotateb * translate4* translate_origin * scale;
    clock_bell.drawSphere(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    int hours, minutes, seconds;
    getCurrentTime(hours, minutes, seconds);
    hours = (hours + 6) % 12;

    
    // second
    shaderActivate(lightingShaderWithTexture);
    model = glm::mat4(1.0f);
    translate_origin = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    scale = glm::scale(model, glm::vec3(0.02, -0.1, 0.01));
    translate = glm::translate(model, glm::vec3(-1.25f, 0.83f + 0.3f + 0.2f, -3.04f));
    glm::mat4 rotateSecond = glm::rotate(model, glm::radians(-(seconds-30)*6.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    translate4 = glm::translate(model, glm::vec3(-1.25f, 0.83f + 0.3f + 0.2f, -3.04f));
    model = alTogether * translate * rotateSecond * translate_origin * scale;
    //model = alTogether * translate * translate_origin * rotateM * translate_origin * scale;
    /*ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);*/
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // minute
    model = glm::mat4(1.0f);
    translate_origin = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    scale = glm::scale(model, glm::vec3(0.02, -0.08, 0.01));
    translate = glm::translate(model, glm::vec3(-1.25f, 0.83f + 0.3f, -3.04f));
    translate2 = glm::translate(model, glm::vec3(0.0f, 0.2f, 0.0f));
    glm::mat4 rotateMinute = glm::rotate(model, glm::radians(-(minutes*60+seconds-30)*0.1f+180.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * translate2 * translate * rotateMinute * translate_origin * scale;
    //model = alTogether * translate * translate_origin * rotateM * translate_origin * scale;
    /*ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);*/
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // hour
    model = glm::mat4(1.0f);
    translate_origin = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    scale = glm::scale(model, glm::vec3(0.02, -0.06, 0.01));
    translate = glm::translate(model, glm::vec3(-1.25f, 0.83f + 0.3f, -3.04f));
    translate2 = glm::translate(model, glm::vec3(0.0f, 0.2f, 0.0f));
    glm::mat4 rotateHour = glm::rotate(model, glm::radians(-(hours*3600+ minutes * 60 + seconds - 30)* (1.0f/120.0f)), glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * translate2 * translate * rotateHour * translate_origin * scale;
    //model = alTogether * translate * translate_origin * rotateM * translate_origin * scale;
    /*ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);*/
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);




    

}


void f_bed_room1_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube)
{
    float baseHeight = 0.01;
    float width = 3.0;
    float length = 3.0;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width + 1, baseHeight, length + 1));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -1.5));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_bed_room1_almari(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& draw_almari)
{
    
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(0.6, 0.0, -4.5));
    scale = glm::scale(model, glm::vec3(0.7, 0.9, 0.3));
    model = alTogether * translate * scale;
    shaderActivate(lightingShaderWithTexture);
    draw_almari.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
    
}


void f_bed_room1_bed(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& bed_sheet, Cube& bed_pillow, Cube& bed_texture, Cube& blanket_texture)
{
    float baseHeight = 0.3;
    float width = 1;
    float length = 2;
    float pillowWidth = 0.3;
    float pillowLength = 0.2;
    float blanketWidth = 0.8;
    float blanketLength = 0.7;
    float headHeight = 0.6;

    //base
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    translate3 = glm::translate(model, glm::vec3(-0.75, 0, -2.7));
    model = alTogether * translate3 * scale * translate;
    shaderActivate(lightingShaderWithTexture);
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    bed_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    //foam
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, baseHeight, 0));
    scale = glm::scale(model, glm::vec3(width, 0.06, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.804, 0.361, 0.361);
    bed_sheet.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    //pillow 1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3((width / 2) - (0.1 + pillowWidth / 2), baseHeight + 1 * 0.06, (length / 2) - (0.025 + pillowWidth / 2)));
    scale = glm::scale(model, glm::vec3(pillowWidth, 0.04, pillowLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 1, 0.647, 0);
    shaderActivate(lightingShaderWithTexture);
    bed_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);

    //pillow 2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3((-width / 2) + (0.1 + pillowWidth / 2), baseHeight + 1 * 0.06, (length / 2) - (0.025 + pillowWidth / 2)));
    scale = glm::scale(model, glm::vec3(pillowWidth, 0.04, pillowLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 1, 0.647, 0);
    bed_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    //blanket
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, baseHeight + 1 * 0.06, -(length / 2 - 0.025) + blanketLength / 2));
    scale = glm::scale(model, glm::vec3(blanketWidth, 0.015, blanketLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.541, 0.169, 0.886);
    blanket_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    //head
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 0, (length / 2 - 0.02 / 2) + 0.02));
    scale = glm::scale(model, glm::vec3(width, headHeight, 0.06));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    bed_texture.drawCubeWithTexture(lightingShaderWithTexture, model);
}


void f_bed_room1_window(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_bed_room1_window)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1, 0.7, 0.05));
    translate = glm::translate(model, glm::vec3(-0.8, 0.5, -4.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    drawing_bed_room1_window.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_bed_room1_bedtable(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_bed_room1_bedtable, Shader& ourShader)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.5, 0.3, 0.2));
    translate = glm::translate(model, glm::vec3(0.0, 0.0, -1.8));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    drawing_bed_room1_bedtable.drawCubeWithTexture(lightingShaderWithTexture, model);

    shaderActivate(ourShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.13, 0.13, 0.13));
    translate = glm::translate(model, glm::vec3(0.17, 0.4, -1.75));
    model = alTogether * translate * scale;
    ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.5f, 0.5f, 0.5f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    
    shaderActivate(lightingShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.03, 0.13, 0.03));
    translate = glm::translate(model, glm::vec3(0.22, 0.3, -1.7));
    model = alTogether * translate * scale;
    drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    //shaderActivate(lightingShaderWithTexture);
    //drawing_bed_room1_bedtable.drawCubeWithTexture(lightingShaderWithTexture, model);

}



// *************************DINING ROOM*************************************
void f_dining_room_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube)
{
    float baseHeight = 1.5;
    float width = 0.01;
    float length = 3.0;

    // right
    float baseHeight1 = baseHeight / 3;
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 rotateYMatrix = glm::rotate(model, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight1, length + 0.1));
    translate = glm::translate(model, glm::vec3(1.6, 2*baseHeight1, -1.5));
    translate2 = glm::translate(model, glm::vec3(3.2, 0.0, -3.1));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);
    
    float baseHeight2 = baseHeight1 * 2;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight2, length + 0.1 - 0.5));
    translate = glm::translate(model, glm::vec3(1.6, 0.0, -1.5));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    
    // left
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // back
    width = 3.0;
    length = 0.01;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width+0.1, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    model = alTogether * translate2 * rotateYMatrix * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
}


void f_dining_room_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube)
{
    float baseHeight = 0.01;
    float width = 3.0;
    float length = 3.0;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length + 1));
    translate = glm::translate(model, glm::vec3(0.65, 0, -1.25));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_dining_table(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& dining_table)
{
    float baseHeight = 0.02;
    float width = 1.2;
    float length = 0.8;

    // Table Top
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0.0, -0.5));
    translate2 = glm::translate(model, glm::vec3(3.25, 0.5, -2.25));
    model = alTogether * translate2 * translate * scale;
    shaderActivate(lightingShaderWithTexture);
    dining_table.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    // Table Legs
    float legWidth = 0.05;
    float legLength = 0.05;

    // Front Left Leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(legWidth, baseHeight - 0.5, legLength));
    translate = glm::translate(model, glm::vec3(-0.5 + 0.1, 0.0, -0.5 + 0.05));
    model = alTogether * translate2 * translate * scale;
    drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);

    // Front Right Leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(legWidth, baseHeight - 0.5, legLength));
    translate = glm::translate(model, glm::vec3(0.5 - legWidth - 0.1 + 0.2, 0.0, -0.5 + 0.05));
    model = alTogether * translate2 * translate * scale;
    drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);

    // Back Left Leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(legWidth, baseHeight - 0.5, legLength));
    translate = glm::translate(model, glm::vec3(-0.5 + 0.1, 0.0, 0.5 - legLength - 0.1 - 0.1));
    model = alTogether * translate2 * translate * scale;
    drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);

    // Back Right Leg
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(legWidth, baseHeight - 0.5, legLength));
    translate = glm::translate(model, glm::vec3(0.5 - legWidth - 0.1 + 0.2, 0.0, 0.5 - legLength - 0.1 - 0.1));
    model = alTogether * translate2 * translate * scale;
    drawCube(cubeVAO, lightingShader, model, 0.24, 0.635, 0.77);
}


void f_dining_frize(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_dining_frize)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.5, 1.0, 0.3));
    translate = glm::translate(model, glm::vec3(1.75, 0.0, -4.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    drawing_dining_frize.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_dining_window(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& kitchen_window)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.02, 0.6, 0.9));
    translate = glm::translate(model, glm::vec3(4.55, 0.6, -2.75));
    model = alTogether * translate * scale;
    kitchen_window.drawCubeWithTexture(lightingShaderWithTexture, model);
}


glm::mat4 transform(float tx, float ty, float tz, float sx, float sy, float sz) {
    glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
    translateMatrix = glm::translate(identityMatrix, glm::vec3(tx, ty, tz));
    rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
    rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
    rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
    scaleMatrix = glm::scale(identityMatrix, glm::vec3(sx, sy, sz));
    model = translateMatrix * scaleMatrix;
    return model;
}


void chair(unsigned int& cubeVAO, Shader& lightingShaderWithTexture, glm::mat4 alTogether, Cube& chair_texture, float tx, float ty, float tz, float rotateY)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 temp = glm::mat4(1.0f);
    glm::mat4 translate_origin = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scaleMatrix = glm::mat4(1.0f);
    scaleMatrix = glm::scale(model, glm::vec3(0.3, 0.3, 0.3));
    translate_origin = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    translate2 = glm::translate(model, glm::vec3(tx + .2, ty - .5, tz + 0.0));
    glm::mat4 rotateYMatrix = glm::rotate(model, glm::radians(rotateY), glm::vec3(0.0f, 1.0f, 0.0f));
    // sit
    temp = transform(tx + .2, ty - .5, tz + 0.0, 1.3, .1, 1);
    model = alTogether * rotateYMatrix * scaleMatrix * temp;
    //model = alTogether * translate2 * rotateYMatrix * translate_origin * temp * scaleMatrix;
    chair_texture.drawCubeWithTexture(lightingShaderWithTexture, model);
    
    // upper
    temp = transform(tx + .2, ty - .5, tz + 1.0, 1.3, 1, .1);
    model = alTogether * rotateYMatrix * scaleMatrix * temp;
    //model = alTogether * translate2 * rotateYMatrix * translate_origin * temp * scaleMatrix;
    chair_texture.drawCubeWithTexture(lightingShaderWithTexture, model);
    
    // lower
    temp = transform(tx + .2, ty - 1.5, tz + 1.0, 1.3, 1, .1);
    model = alTogether * rotateYMatrix * scaleMatrix * temp;
    //model = alTogether * translate2 * rotateYMatrix * translate_origin * temp * scaleMatrix;
    chair_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    // lower
    temp = transform(tx + .2, ty - 1.5, tz + 0.0, 1.3, 1, .1);
    model = alTogether * rotateYMatrix * scaleMatrix * temp;
    //model = alTogether * translate2 * rotateYMatrix * translate_origin * temp * scaleMatrix;
    chair_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

}


void f_dining_light(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Shader& ourShader)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // point light
    shaderActivate(ourShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.15, 0.15));
    translate = glm::translate(model, glm::vec3(1.6f, 1.3f, -3.1f));
    model = alTogether * translate * scale;
    ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    shaderActivate(lightingShader);
}
// **********************************       Kitchen           ***********************************
void f_kitchen_surface(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& kitchen_surface, Cube& kitchen_surface_top, Cube& kitchen_cupboard1, Cube& kitchen_cupboard2, Cube& kitchen_back_texture)
{
    // ground
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, 0.5, 0.4));
    translate = glm::translate(model, glm::vec3(2.75, 0.0, -4.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    kitchen_surface.drawCubeWithTexture(lightingShaderWithTexture, model);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, 0.02, 0.4));
    translate = glm::translate(model, glm::vec3(2.75, 0.5, -4.5));
    model = alTogether * translate * scale;
    kitchen_surface_top.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.4, 0.5, 0.9));
    translate = glm::translate(model, glm::vec3(4.25, 0.0, -4.1));
    model = alTogether * translate * scale;
    kitchen_surface.drawCubeWithTexture(lightingShaderWithTexture, model);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.4, 0.02, 0.9));
    translate = glm::translate(model, glm::vec3(4.25, 0.5, -4.1));
    model = alTogether * translate * scale;
    kitchen_surface_top.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, 0.5, 0.4));
    translate = glm::translate(model, glm::vec3(2.25, 0.0, -3.6));
    model = alTogether * translate * scale;
    kitchen_surface.drawCubeWithTexture(lightingShaderWithTexture, model);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, 0.02, 0.4));
    translate = glm::translate(model, glm::vec3(2.25, 0.5, -3.6));
    model = alTogether * translate * scale;
    kitchen_surface_top.drawCubeWithTexture(lightingShaderWithTexture, model);

    // cupboard
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, 0.5, 0.4));
    translate = glm::translate(model, glm::vec3(2.75, 1.0, -4.5));
    model = alTogether * translate * scale;
    kitchen_cupboard1.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.4, 0.5, 0.9));
    translate = glm::translate(model, glm::vec3(4.25, 1.0, -4.1));
    model = alTogether * translate * scale;
    kitchen_cupboard2.drawCubeWithTexture(lightingShaderWithTexture, model);

    // kitchen_back
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1.5, 0.4, 0.02));
    translate = glm::translate(model, glm::vec3(2.75, 0.5, -4.5));
    model = alTogether * translate * scale;
    kitchen_back_texture.drawCubeWithTexture(lightingShaderWithTexture, model);


    shaderActivate(lightingShader);
}


void f_kitchen_stove(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& kitchen_surface)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.6, 0.06, 0.4));
    translate = glm::translate(model, glm::vec3(3.3, 0.52, -3.6));
    model = alTogether * translate * scale;
    kitchen_surface.drawCubeWithTexture(lightingShaderWithTexture, model);
    
}
// **********************************       BED ROOM 2        ************************************
void f_bed_room2_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube)
{
    float baseHeight = 1.5;
    float width = 0.01;
    float length = 3.0;

    // right
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(1.6, 0.0, -1.5));
    translate2 = glm::translate(model, glm::vec3(3.0, 0.0, 0.0));
    model = alTogether * translate2 * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);


    // back
    // back top
    float baseHeight1 = baseHeight / 3;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(length + 0.1, baseHeight1, width + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 2 * baseHeight1, 1.5));
    glm::mat4 translate3 = glm::translate(model, glm::vec3(0.0, 0.0, -1.75));
    model = alTogether * translate3 * translate2 * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    // front bottom
    float baseHeight2 = baseHeight1 * 2;
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(length - 0.5, baseHeight2, width + 0.1));
    translate = glm::translate(model, glm::vec3(-0.9, 0.0, 1.5));
    model = alTogether * translate3 * translate2 * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    // back
    width = 3.0;
    length = 0.01;
    //model = glm::mat4(1.0f);
    //scale = glm::scale(model, glm::vec3(width, baseHeight, length + 0.1));
    //translate = glm::translate(model, glm::vec3(-1.5, 0.0, -1.5));
    //model = alTogether * translate * scale;
    ////drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //

    // front
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width + 0.1, baseHeight, length + 0.1));
    translate = glm::translate(model, glm::vec3(-1.5, 0.0, 1.5));
    model = alTogether * translate2 * translate * scale;
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    

    shaderActivate(lightingShader);

}


void f_bed_room2_floor(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& cube)
{
    float baseHeight = 0.01;
    float width = 3.0;
    float length = 3.0;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length + 1));
    translate = glm::translate(model, glm::vec3(0.65, 0, -0.5));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_bed_room2_book(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_dining_frize)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.8, 0.6, 0.2));
    translate = glm::translate(model, glm::vec3(1.75, 0.5, 1.3));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    shaderActivate(lightingShaderWithTexture);
    drawing_dining_frize.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void f_bed_room2_light(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Shader& ourShader)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // point light
    shaderActivate(ourShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.15, 0.15));
    translate = glm::translate(model, glm::vec3(1.6f, 1.3f, 0.5f));
    model = alTogether * translate * scale;
    ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    shaderActivate(lightingShader);
}


void f_bed_room2_bed(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& bed_sheet, Cube& bed_pillow, Cube& bed_texture, Cube& blanket_texture)
{
    float baseHeight = 0.3;
    float width = 0.7;
    float length = 2;
    float pillowWidth = 0.25;
    float pillowLength = 0.2;
    float blanketWidth = 0.7;
    float blanketLength = 0.5;
    float headHeight = 0.6;

    //base
    glm::mat4 translate_origin = glm::mat4(1.0f);
    glm::mat4 translate_location = glm::mat4(1.0f);
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 rotateM = glm::rotate(model, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    translate_origin = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    translate_location = glm::translate(model, glm::vec3(6.3, 0.0, -0.45));
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    translate3 = glm::translate(model, glm::vec3(-0.75, 0, -2.7));
    model = alTogether * translate_location * rotateM * translate3 * scale * translate;
    shaderActivate(lightingShaderWithTexture);
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    bed_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    //foam
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, baseHeight, 0));
    scale = glm::scale(model, glm::vec3(width, 0.06, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate_location * rotateM * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.804, 0.361, 0.361);
    bed_sheet.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    //pillow 1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3((width / 2) - (0.1 + pillowWidth / 2), baseHeight + 1 * 0.06, (length / 2) - (0.025 + pillowWidth / 2)));
    scale = glm::scale(model, glm::vec3(pillowWidth, 0.04, pillowLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate_location * rotateM * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 1, 0.647, 0);
    shaderActivate(lightingShaderWithTexture);
    bed_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);

    //pillow 2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3((-width / 2) + (0.1 + pillowWidth / 2), baseHeight + 1 * 0.06, (length / 2) - (0.025 + pillowWidth / 2)));
    scale = glm::scale(model, glm::vec3(pillowWidth, 0.04, pillowLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate_location * rotateM * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 1, 0.647, 0);
    bed_pillow.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);

    //blanket
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, baseHeight + 1 * 0.06, -(length / 2 - 0.025) + blanketLength / 2));
    scale = glm::scale(model, glm::vec3(blanketWidth, 0.015, blanketLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate_location * rotateM * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.541, 0.169, 0.886);
    blanket_texture.drawCubeWithTexture(lightingShaderWithTexture, model);

    //head
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 0, (length / 2 - 0.02 / 2) + 0.02));
    scale = glm::scale(model, glm::vec3(width, headHeight, 0.06));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate_location * rotateM * translate3 * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    bed_texture.drawCubeWithTexture(lightingShaderWithTexture, model);
}


void f_bed_room2_dressing_table(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& dressing_right, Cube& dressing_bottom, Cube& dressing_mirror)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // right part
    shaderActivate(lightingShaderWithTexture);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.2, 1.0, 0.3));
    translate = glm::translate(model, glm::vec3(4.4f, 0.0f, 1.2f));
    model = alTogether * translate * scale;
    dressing_right.drawCubeWithTexture(lightingShaderWithTexture, model);

    // bottom part
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.2, 0.3, 0.3));
    translate = glm::translate(model, glm::vec3(4.4f, 0.0f, 0.9f));
    model = alTogether * translate * scale;
    dressing_bottom.drawCubeWithTexture(lightingShaderWithTexture, model);

    // mirror
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.6, 0.3));
    translate = glm::translate(model, glm::vec3(4.55f, 0.35f, 0.9f));
    model = alTogether * translate * scale;
    dressing_mirror.drawCubeWithTexture(lightingShaderWithTexture, model);


}
// *********************************   ROOF    *****************************************************
void f_roof(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& drawing_roof)
{
    float baseHeight = 0.01;
    float width = 3.0;
    float length = 3.0;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2*(width + 1), baseHeight, 2*(length + 1)));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    translate2 = glm::translate(model, glm::vec3(1.5, 1.5, -1.5));
    model = alTogether * translate2 * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    drawing_roof.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}

// *************************BATHROOM****************************
void f_bathroom_wall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& bathroom_top, Cube& bathroom_door, Cube& bathroom_tiles)
{

    float baseHeight = 1.5;
    float width = 0.01;
    float length = 3.0;
    float base_top = baseHeight / 3.0;
    shaderActivate(lightingShaderWithTexture);

    // top
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, base_top, length + 0.1 - 1.75));
    translate = glm::translate(model, glm::vec3(2.5, 2*base_top, -1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    bathroom_top.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left_door
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, 2*base_top, (length + 0.1 - 1.75)/2.0));
    translate = glm::translate(model, glm::vec3(2.5, 0.0, -1.5 + bathroom_door_translate));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_door.drawCubeWithTexture(lightingShaderWithTexture, model);

    
    // right_door
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, 2 * base_top, (length + 0.1 - 1.75) / 2.0));
    translate = glm::translate(model, glm::vec3(2.5, 0.0, -1.5+ (length + 0.1 - 1.75) / 2.0));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_door.drawCubeWithTexture(lightingShaderWithTexture, model);

    // back_wall
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, 2 * base_top, length + 0.1 - 1.75));
    translate = glm::translate(model, glm::vec3(2.5 + 2.0, 0.0, -1.5));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_tiles.drawCubeWithTexture(lightingShaderWithTexture, model);

    // left_wall
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(length + 0.1 - 1.25, 2 * base_top, width));
    translate = glm::translate(model, glm::vec3(2.6, 0.0, -1.5+0.15));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_tiles.drawCubeWithTexture(lightingShaderWithTexture, model);

    // right_wall
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(length + 0.1 - 1.25, 2 * base_top, width));
    translate = glm::translate(model, glm::vec3(2.6, 0.0, -0.3));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_tiles.drawCubeWithTexture(lightingShaderWithTexture, model);

    // top_wall
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, width, 2*max_bathroom_door_translate));
    translate = glm::translate(model, glm::vec3(2.6, 2*base_top, -1.65 + 0.15));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_tiles.drawCubeWithTexture(lightingShaderWithTexture, model);

    // bottom wall
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(2.0, width, 2 * max_bathroom_door_translate));
    translate = glm::translate(model, glm::vec3(2.6, 0.1, -1.65 + 0.15));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.3, 0.3, 0.3);
    shaderActivate(lightingShaderWithTexture);
    bathroom_tiles.drawCubeWithTexture(lightingShaderWithTexture, model);


    shaderActivate(lightingShader);
}

void f_bathroom_light(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Shader& ourShader)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // point light
    shaderActivate(ourShader);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.05, 0.15, 0.15));
    translate = glm::translate(model, glm::vec3(2.5 + 1.9, 0.8, -0.9));
    model = alTogether * translate * scale;
    ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    shaderActivate(lightingShader);
}

void f_toilet(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& toilet)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // toilet
    shaderActivate(lightingShaderWithTexture);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.4, 0.05, 0.6));
    translate = glm::translate(model, glm::vec3(2.5 + 1.5, 0.1, -0.95));
    model = alTogether * translate * scale;
    /*ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);*/
    toilet.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}

void f_bathroom_bucket(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, Cube& water)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);

    // toilet
    shaderActivate(lightingShaderWithTexture);
    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.25, 0.25, 0.5));
    translate = glm::translate(model, glm::vec3(2.5 + 1.0, 0.05, -1.1));
    translate2 = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    glm::mat4 rotateM = glm::rotate(model, glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));

    model = alTogether * translate * rotateM * translate2 * scale;
    /*ourShader.setMat4("model", model);
    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);*/
    //toilet.drawCubeWithTexture(lightingShaderWithTexture, model);
    bucket->draw(lightingShaderWithTexture, model);
    //shaderActivate(lightingShader);

    model = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.3, 0.21, 0.5));
    translate = glm::translate(model, glm::vec3(2.5 + 0.75, 0.3, -0.95));
    model = alTogether * translate * rotateM * translate2 * scale;
    water.drawCubeWithTexture(lightingShaderWithTexture, model);
    shaderActivate(lightingShader);
}


void getCurrentTime(int& hours, int& minutes, int& seconds) {
    time_t currentTime = time(nullptr); // Get current UNIX timestamp
    struct tm* timeinfo;
    timeinfo = localtime(&currentTime);

    seconds = timeinfo->tm_sec;
    minutes = timeinfo->tm_min;
    hours = timeinfo->tm_hour;
}
// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        camera.ProcessKeyboard(FORWARD, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        camera.ProcessKeyboard(LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        camera.ProcessKeyboard(RIGHT, deltaTime);
    }

    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
        camera.ProcessKeyboard(UP, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
        camera.ProcessKeyboard(DOWN, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS) {
        camera.ProcessKeyboard(P_UP, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) {
        camera.ProcessKeyboard(P_DOWN, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) {
        camera.ProcessKeyboard(Y_LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS) {
        camera.ProcessKeyboard(Y_RIGHT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) {
        camera.ProcessKeyboard(R_LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_J) == GLFW_PRESS) {
        camera.ProcessKeyboard(R_RIGHT, deltaTime);
    }

    /*if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
    {
        if (rotateAxis_X) rotateAngle_X -= 0.1;
        else if (rotateAxis_Y) rotateAngle_Y -= 0.1;
        else rotateAngle_Z -= 0.1;
    }*/
    /*if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) translate_Y += 0.001;*/
    //if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS) translate_Y -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) translate_X += 0.001;
    /*if (glfwGetKey(window, GLFW_KEY_J) == GLFW_PRESS) translate_X -= 0.001;*/
    if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) translate_Z += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) translate_Z -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS) scale_X += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS) scale_X -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS) scale_Y += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS) scale_Y -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS) scale_Z += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) scale_Z -= 0.001;

    if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
    {
        rotateAngle_X += 0.1;
        rotateAxis_X = 1.0;
        rotateAxis_Y = 0.0;
        rotateAxis_Z = 0.0;
    }
    if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS)
    {
        rotateAngle_Y += 0.1;
        rotateAxis_X = 0.0;
        rotateAxis_Y = 1.0;
        rotateAxis_Z = 0.0;
    }
    if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
    {
        rotateAngle_Z += 0.1;
        rotateAxis_X = 0.0;
        rotateAxis_Y = 0.0;
        rotateAxis_Z = 1.0;
    }

    /*if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS)
    {
        eyeX += 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }*/
    /*if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS)
    {
        eyeX -= 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }*/
    /*if (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS)
    {
        eyeZ += 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }*/
    if (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS)
    {
        eyeZ -= 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
    {
        eyeY += 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    /*if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
    {
        eyeY -= 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }*/
    // front door
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    {
        f_door += 1;
        f_door = min(70.0f, f_door);
    }
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
    {
        f_door -= 1;
        f_door = max(0.0f, f_door);
    }
    // back door
    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
    {
        b1_door += 1;
        b1_door = min(80.0f, b1_door);
    }
    if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
    {
        b1_door -= 1;
        b1_door = max(0.0f, b1_door);
    }
    // bathroom door
    if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
    {
        bathroom_door_translate += 0.01;
        bathroom_door_translate = min(max_bathroom_door_translate, bathroom_door_translate);
    }
    if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
    {
        bathroom_door_translate -= 0.01;
        bathroom_door_translate = max(0.0f, bathroom_door_translate);
    }


    // fan
    if (glfwGetKey(window, GLFW_KEY_0) == GLFW_PRESS)
    {
        rotateFan += 5.0f;
    }
    
    if (glfwGetKey(window, GLFW_KEY_7) == GLFW_PRESS)
    {
        pointlight1.turnAmbientOn();
        pointlight3.turnAmbientOn();
        drawing_light.turnAmbientOn();
        bed_room1_light.turnAmbientOn();
        dining_light.turnAmbientOn();
        bed_room2_light.turnAmbientOn();
        bathroom_light.turnAmbientOn();

        pointlight1.turnDiffuseOff();
        pointlight3.turnDiffuseOff();
        drawing_light.turnDiffuseOff();
        bed_room1_light.turnDiffuseOff();
        dining_light.turnDiffuseOff();
        bed_room2_light.turnDiffuseOff();
        bathroom_light.turnDiffuseOff();

        pointlight1.turnSpecularOff();
        pointlight3.turnSpecularOff();
        drawing_light.turnSpecularOff();
        bed_room1_light.turnSpecularOff();
        dining_light.turnSpecularOff();
        bed_room2_light.turnSpecularOff();
        bathroom_light.turnSpecularOff();
    }
    
    if (glfwGetKey(window, GLFW_KEY_8) == GLFW_PRESS)
    {
        pointlight1.turnDiffuseOn();
        pointlight3.turnDiffuseOn();
        drawing_light.turnDiffuseOn();
        bed_room1_light.turnDiffuseOn();
        dining_light.turnDiffuseOn();
        bed_room2_light.turnDiffuseOn();
        bathroom_light.turnDiffuseOn();

        pointlight1.turnAmbientOff();
        pointlight3.turnAmbientOff();
        drawing_light.turnAmbientOff();
        bed_room1_light.turnAmbientOff();
        dining_light.turnAmbientOff();
        bed_room2_light.turnAmbientOff();
        bathroom_light.turnAmbientOff();

        pointlight1.turnSpecularOff();
        pointlight3.turnSpecularOff();
        drawing_light.turnSpecularOff();
        bed_room1_light.turnSpecularOff();
        dining_light.turnSpecularOff();
        bed_room2_light.turnSpecularOff();
        bathroom_light.turnSpecularOff();
    }

    if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS)
    {
        pointlight1.turnDiffuseOff();
        pointlight3.turnDiffuseOff();
        drawing_light.turnDiffuseOff();
        bed_room1_light.turnDiffuseOff();
        dining_light.turnDiffuseOff();
        bed_room2_light.turnDiffuseOff();
        bathroom_light.turnDiffuseOff();

        pointlight1.turnAmbientOff();
        pointlight3.turnAmbientOff();
        drawing_light.turnAmbientOff();
        bed_room1_light.turnAmbientOff();
        dining_light.turnAmbientOff();
        bed_room2_light.turnAmbientOff();
        bathroom_light.turnAmbientOff();

        pointlight1.turnSpecularOn();
        pointlight3.turnSpecularOn();
        drawing_light.turnSpecularOn();
        bed_room1_light.turnSpecularOn();
        dining_light.turnSpecularOn();
        bed_room2_light.turnSpecularOn();
        bathroom_light.turnSpecularOn();
    }

    if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS)
    {
        pointlight1.turnDiffuseOn();
        pointlight3.turnDiffuseOn();
        drawing_light.turnDiffuseOn();
        bed_room1_light.turnDiffuseOn();
        dining_light.turnDiffuseOn();
        bed_room2_light.turnDiffuseOn();
        bathroom_light.turnDiffuseOn();

        pointlight1.turnAmbientOn();
        pointlight3.turnAmbientOn();
        drawing_light.turnAmbientOn();
        bed_room1_light.turnAmbientOn();
        dining_light.turnAmbientOn();
        bed_room2_light.turnAmbientOn();
        bathroom_light.turnAmbientOn();

        pointlight1.turnSpecularOn();
        pointlight3.turnSpecularOn();
        drawing_light.turnSpecularOn();
        bed_room1_light.turnSpecularOn();
        dining_light.turnSpecularOn();
        bed_room2_light.turnSpecularOn();
        bathroom_light.turnSpecularOn();
    }


   // if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
   // {
   //     if (onOffPointToggle)
   //     {
   //         pointlight1.turnOff();
   //         
   //         onOffPointToggle = false;
   //     }
   //     else
   //     {
   //         pointlight1.turnOn();
   //       
   //         onOffPointToggle = true;
   //     }
   //    // pointlight3.turnOff();
   //    // pointlight4.turnOff();

   // }
   // 

   // if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
   // {
   //     
   //     if (onOffSpotToggle)
   //     {
   //        
   //         pointlight2.turnOff();
   //         onOffSpotToggle = false;
   //     }
   //     else
   //     {
   //         pointlight2.turnOn();
   //         onOffSpotToggle = true;
   //     }
   // }

   // if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
   // {

   //     if (onOffDirectToggle)
   //     {

   //         pointlight3.turnOff();
   //         onOffDirectToggle = false;
   //     }
   //     else
   //     {
   //         pointlight3.turnOn();
   //         onOffDirectToggle = true;
   //     }
   // }
   // 
   // if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
   // {
   //     pointlight1.turnAmbientOn();
   //     pointlight2.turnAmbientOn();
   //    // pointlight3.turnAmbientOn();
   //    // pointlight4.turnAmbientOn();
   // }
   // if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
   // {
   //     pointlight1.turnAmbientOff();
   //     pointlight2.turnAmbientOff();
   //   //  pointlight3.turnAmbientOff();
   //   //  pointlight4.turnAmbientOff();
   // }
   // if (glfwGetKey(window, GLFW_KEY_7) == GLFW_PRESS)
   // {
   //     pointlight1.turnDiffuseOn();
   //     pointlight2.turnDiffuseOn();
   //  //   pointlight3.turnDiffuseOn();
   // //    pointlight4.turnDiffuseOn();
   // }
   // if (glfwGetKey(window, GLFW_KEY_8) == GLFW_PRESS)
   // {
   //     pointlight1.turnDiffuseOff();
   //     pointlight2.turnDiffuseOff();
   ////     pointlight3.turnDiffuseOff();
   // //    pointlight4.turnDiffuseOff();
   // }
   // if (glfwGetKey(window, GLFW_KEY_9) == GLFW_PRESS)
   // {
   //     pointlight1.turnSpecularOn();
   //     pointlight2.turnSpecularOn();
   // //    pointlight3.turnSpecularOn();
   // //    pointlight4.turnSpecularOn();
   // }
   // if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
   // {
   //     pointlight1.turnSpecularOff();
   //     pointlight2.turnSpecularOff();
   ////     pointlight3.turnSpecularOff();
   // //    pointlight4.turnDiffuseOff();
   // }
    if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS)
    {
        //pointlight2.turnOn();
        drawing_light.turnOn();
        bed_room1_light.turnOn();
        dining_light.turnOn();
        bed_room2_light.turnOn();
        bathroom_light.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS)
    {
        //pointlight2.turnOff();
        drawing_light.turnOff();
        bed_room1_light.turnOff();
        dining_light.turnOff();
        bed_room2_light.turnOff();
        bathroom_light.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS)
    {
        pointlight3.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS)
    {
        pointlight3.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS)
    {
        pointlight1.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
    {
        pointlight1.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    //if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    //{
    //    /*if (diffuseToggle)
    //    {*/
    //    if (pointlight1.isOn())
    //        pointlight1.turnAmbientOn();
    //    if (pointlight2.isOn())
    //        pointlight2.turnAmbientOn();
    //    if (pointlight3.isOn())
    //        pointlight3.turnAmbientOn();
    //    //pointlight4.turnDiffuseOn();
    //    //diffuseToggle = !diffuseToggle;
    ////}
    //}
    //if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
    //{
    //    /*if (diffuseToggle)
    //    {*/
    //    if (pointlight1.isOn())
    //        pointlight1.turnAmbientOff();
    //    if (pointlight2.isOn())
    //        pointlight2.turnAmbientOff();
    //    if (pointlight3.isOn())
    //        pointlight3.turnAmbientOff();
    //    //pointlight4.turnDiffuseOff();
    //    //diffuseToggle = !diffuseToggle;
    ////}
    //}

    //if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
    //{
    //    /*if (diffuseToggle)
    //    {*/
    //    if (pointlight1.isOn())
    //        pointlight1.turnDiffuseOn();
    //    if (pointlight2.isOn())
    //        pointlight2.turnDiffuseOn();
    //    if (pointlight3.isOn())
    //        pointlight3.turnDiffuseOn();
    //    //pointlight4.turnAmbientOn();
    //    //diffuseToggle = !diffuseToggle;
    //    //}
    //}
    //if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
    //{
    //    /*if (diffuseToggle)
    //    {*/
    //    if (pointlight1.isOn())
    //        pointlight1.turnDiffuseOff();
    //    if (pointlight2.isOn())
    //        pointlight2.turnDiffuseOff();
    //    if (pointlight3.isOn())
    //        pointlight3.turnDiffuseOff();
    //    //diffuseToggle = !diffuseToggle;
    //    //}
    //}


    if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnSpecularOn();
        if (pointlight2.isOn())
            pointlight2.turnSpecularOn();
        if (pointlight3.isOn())
            pointlight3.turnSpecularOn();
        //pointlight4.turnSpecularOn();
        //diffuseToggle = !diffuseToggle;
        //}
    }
    if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        /*cout << "1 " << pointlight1.isOn() << endl;
        cout << pointlight2.isOn() << endl;
        cout << pointlight3.isOn() << endl;*/
        if (pointlight1.isOn())
            pointlight1.turnSpecularOff();
        if (pointlight2.isOn())
            pointlight2.turnSpecularOff();
        if (pointlight3.isOn())
            pointlight3.turnSpecularOff();
        //pointlight4.turnSpecularOff();
        //diffuseToggle = !diffuseToggle;
        //}
    }
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}


// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(static_cast<float>(yoffset));
}

unsigned int loadTexture(char const* path, GLenum textureWrappingModeS, GLenum textureWrappingModeT, GLenum textureFilteringModeMin, GLenum textureFilteringModeMax)
{
    unsigned int textureID;
    glGenTextures(1, &textureID);

    int width, height, nrComponents;
    stbi_set_flip_vertically_on_load(true);
    unsigned char* data = stbi_load(path, &width, &height, &nrComponents, 0);
    if (data)
    {
        GLenum format;
        if (nrComponents == 1)
            format = GL_RED;
        else if (nrComponents == 3)
            format = GL_RGB;
        else if (nrComponents == 4)
            format = GL_RGBA;

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, textureWrappingModeS);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, textureWrappingModeT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, textureFilteringModeMin);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, textureFilteringModeMax);

        stbi_image_free(data);
    }
    else
    {
        std::cout << "Texture failed to load at path: " << path << std::endl;
        stbi_image_free(data);
    }

    return textureID;
}
