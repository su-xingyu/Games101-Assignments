#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.

        float dx = (end.x-start.x) / (num_nodes-1);
        float dy = (end.y-start.y) / (num_nodes-1);

        float x = start.x;
        float y = start.y;

        for (int i = 0; i < num_nodes; i++) {

            Mass* m_ptr = new Mass(Vector2D(x, y), node_mass, false);
            masses.push_back(m_ptr);

            if (i > 0) {
                Spring* s_ptr = new Spring(masses[i-1], masses[i], k);
                springs.push_back(s_ptr);
            }

            x += dx;
            y += dy;
        }

//        Comment-in this part when you implement the constructor
        for (auto &i : pinned_nodes) {
            masses[i]->pinned = true;
        }
        
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            Vector2D l = s->m2->position - s->m1->position;
            Vector2D f_a2b = s->k * (l.norm()-s->rest_length) * l.unit();

            if (!s->m1->pinned) {
                s->m1->forces += f_a2b;
            }
            if (!s->m2->pinned) {
                s->m2->forces += -f_a2b;
            }
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                m->forces += m->mass * gravity;
                m->velocity += m->forces / m->mass * delta_t;
                m->position += m->velocity * delta_t;
                // TODO (Part 2): Add global damping
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            Vector2D l = s->m2->position - s->m1->position;
            Vector2D f_a2b = s->k * (l.norm()-s->rest_length) * l.unit();

            if (!s->m1->pinned) {
                s->m1->forces += f_a2b;
            }
            if (!s->m2->pinned) {
                s->m2->forces += -f_a2b;
            }
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                Vector2D temp_position = m->position;
                // TODO (Part 3.1): Set the new position of the rope mass
                m->forces += m->mass * gravity;
                Vector2D a = m->forces / m->mass;

                // TODO (Part 4): Add global Verlet damping
                float damping_factor = 0.00005;
                m->position += (1-damping_factor) * (m->position-m->last_position) + a * delta_t * delta_t;
                m->last_position = temp_position;
            }
            m->forces = Vector2D(0, 0);
        }
    }
}
