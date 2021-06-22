/*
 *       /\        Matteo Cicuttin (C) 2016, 2017, 2018
 *      /__\       matteo.cicuttin@enpc.fr
 *     /_\/_\      École Nationale des Ponts et Chaussées - CERMICS
 *    /\    /\
 *   /__\  /__\    DISK++, a template library for DIscontinuous SKeletal
 *  /_\/_\/_\/_\   methods.
 *
 * This file is copyright of the following authors:
 * Nicolas Pignet  (C) 2019                     nicolas.pignet@enpc.fr
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * If you use this code or parts of it for scientific publications, you
 * are required to cite it as following:
 *
 * Hybrid High-Order methods for finite elastoplastic deformations
 * within a logarithmic strain framework.
 * M. Abbas, A. Ern, N. Pignet.
 * International Journal of Numerical Methods in Engineering (2019)
 * 120(3), 303-327
 * DOI: 10.1002/nme.6137
 */

#pragma once

#include <iostream>
#include <list>
#include <vector>

namespace disk
{

namespace mechanics
{

/**
 * @brief Represent a time step on the time interval [start_time, end_time]
 *
 * @tparam T scalar type
 */
template<typename T>
class TimeStep
{
  private:
    T      m_start_time, m_end_time;
    size_t m_level;

  public:
    TimeStep() : m_start_time(T(0)), m_end_time(T(0)), m_level(0) {}

    /**
     * @brief Construct a new TimeStep object
     *
     * @param start_time beginning of the time step
     * @param end_time end of the time step
     * @param level level of time refining
     */
    TimeStep(const T start_time, const T end_time, const size_t level) :
      m_start_time(start_time), m_end_time(end_time), m_level(level)
    {
    }

    /**
     * @brief Return the start of the time step
     *
     * @return T  start time
     */
    T
    start_time(void) const
    {
        return m_start_time;
    }

    /**
     * @brief Return the end of the time step
     *
     * @return T end time
     */
    T
    end_time(void) const
    {
        return m_end_time;
    }

    /**
     * @brief Return the level of time refining
     *
     * @return size_t level
     */
    size_t
    level(void) const
    {
        return m_level;
    }
};

/**
 * @brief This class ListOfTimeStep allows to manage the different time step for the simulation
 *
 * @tparam T scalar type
 */
template<typename T>
class ListOfTimeStep
{
  private:
    std::list<TimeStep<T>> list_steps;
    size_t                 n_time_step_comp;
    T                      user_end_time;

    template<typename I>
    void
    addTimeStepping(const T end_time, const I n_step, const size_t level)
    {
        T start_time = T(0);

        if (!list_steps.empty())
        {
            start_time = list_steps.back().end_time();
        }

        const T delta_t = (end_time - start_time) / T(n_step);

        for (I i = 0; i < n_step; i++)
        {
            if ((start_time + (i + 1) * delta_t) <= user_end_time)
            {
                const TimeStep<T> step(start_time + i * delta_t, start_time + (i + 1) * delta_t, level);
                list_steps.push_back(step);
            }
        }
    }

  public:
    ListOfTimeStep() : n_time_step_comp(0), user_end_time(0)
    {
        list_steps.clear();
    }

    ListOfTimeStep(const T end_time, const size_t n_step) : n_time_step_comp(0), user_end_time(end_time)
    {
        this->addTimeStepping(end_time, n_step, 0);
    }

    template<typename I>
    ListOfTimeStep(const std::vector<std::pair<T, I>> list_of_time_step) : n_time_step_comp(0)
    {
        user_end_time = list_of_time_step[list_of_time_step.size()-1].first;
        for (auto& [time, n_step] : list_of_time_step)
        {
            this->addTimeStepping(time, n_step, 0);
        }
    }

    template<typename I>
    ListOfTimeStep(const std::vector<std::pair<T, I>> list_of_time_step, const T final_time) : n_time_step_comp(0)
    {
        user_end_time = final_time;
        for (auto& [time, n_step] : list_of_time_step) { this->addTimeStepping(time, n_step, 0); }
    }

    /**
     * @brief Return a boolean to known if there is still some time step
     *
     * @return true No more time step
     */
    bool
    empty(void) const
    {
        return list_steps.empty();
    }

    /**
     * @brief Return the number of time step
     *
     * @return size_t number of time step
     */
    size_t
    numberOfTimeStep(void) const
    {
        return list_steps.size() + n_time_step_comp;
    }

    /**
     * @brief Return the number of time step
     *
     * @return size_t number of time step
     */
    size_t
    numberOfRemainingTimeStep(void) const
    {
        return list_steps.size();
    }

    /**
     * @brief Return the number of time step already realized
     *
     * @return size_t number of time step realisized
     */
    size_t
    numberOfTimeStepRealized(void) const
    {
        return n_time_step_comp;
    }

    /**
     * @brief Return the index (n-th)-time step
     *
     * @return size_t index of current time step
     */
    size_t
    indexOfCurrentTimeStep(void) const
    {
        return n_time_step_comp + 1;
    }

    /**
     * @brief Get the Current Time Step object
     *
     * @return TimeStep<T> Current time step
     */
    TimeStep<T>
    getCurrentTimeStep(void) const
    {
        return list_steps.front();
    }

    /**
     * @brief Remove the current time step
     *
     */
    void
    removeCurrentTimeStep(void)
    {
        list_steps.pop_front();
        n_time_step_comp++;
    }

    /**
     * @brief Split the Current time step in two equal sub-time step. This allows to refine the time step
     *
     */
    void
    splitCurrentTimeStep(void)
    {
        const TimeStep<T> CurrentStep = this->getCurrentTimeStep();

        const T start_time = CurrentStep.start_time();
        const T end_time   = CurrentStep.end_time();

        const T new_time = start_time + (end_time - start_time) / T(2);

        const TimeStep<T> newStep1(start_time, new_time, CurrentStep.level() + 1);
        const TimeStep<T> newStep2(new_time, end_time, CurrentStep.level() + 1);

        list_steps.pop_front();
        list_steps.push_front(newStep2);
        list_steps.push_front(newStep1);
    }

    /**
     * @brief Print informations about current time step
     *
     */
    void
    printCurrentTimeStep(void) const
    {
        const TimeStep<T> current_step = this->getCurrentTimeStep();
        const T           current_time = current_step.end_time();

        std::cout << "------------------------------------------------------------------------"
                     "----------------------"
                  << std::endl;
        std::cout << "************************ Time : " << current_time
                  << " sec (step: " << this->indexOfCurrentTimeStep() << "/" << this->numberOfTimeStep()
                  << ", sublevel: " << current_step.level() << " ) *************************|" << std::endl;
    }
};
}
}