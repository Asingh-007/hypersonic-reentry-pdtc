#ifndef CONTROL_INPUTS_H
#define CONTROL_INPUTS_H

struct ThrustVectorControlInputs {
    float thrust = 0.0f; // thrust magnitude in Newtons
    float pitch = 0.0f; // pitch angle in radians
    float yaw = 0.0f; // yaw angle in radians

    ThrustVectorControlInputs(float thrust, float pitch, float yaw)
        : thrust(thrust), pitch(pitch), yaw(yaw) {}
};

#endif // CONTROL_INPUTS_H