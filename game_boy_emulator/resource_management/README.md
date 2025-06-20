# Resource Management

Resource management is a very important topic in software
development. Here are a few examples of how I managed it in my code.

(I'll admit there are still memory leaks in my software, but I'm quite
certain they are purely from the SDL2 audio/video library. Even
official examples from the SDL2 library have similar leaks, so I'm not
going to dwell on those.)

## Battery

The battery is a save file that, during gameplay, can get written to
repeatedly. My implementation of [battery](battery.hpp) follows the
RAII paradigm. It opens the save file in a constructor so it is ready
to go for writing. When a battery object gets destroyed, the file
closes automatically.  Given how frequently data can be written during
gameplay, keeping this file open the whole time, rather than
repeatedly opening and closing it, can help reduce the computational
cost of saving data.

## Instructions

The Game Boy CPU includes 500 instructions. My implementation required
storing them in some sort of vector. I did this by using arrays which
get freed automatically when emulator closes. Additionally, I made the
instructions unique pointers to ensure they get freed as well. The
instructions get returned through the ```get_instruction()``` function
as references to ensure there is no confusion about ownership.

## Global Variables

The emulator includes 3 processing units and various other objects that
interact. I decided to make them all global, by constructing them in a
namespace (see [global.cpp](global.cpp)). Then in
[main.cpp](main.cpp), after they are all guaranteed to be constructed,
they are initialized. These objects are all guaranteed to be unique,
their global access reduces passing them as additional arguments,
sharing them as pointers, etc., and when the emulator closes they are
all guaranteed to be destructed.

A huge red flag with my design is how all the objects are in a single
file, which can have a significant impact on compilation time! Any
future work on this project will likely start by reorganizing this
extensively to accelerate compilation.
