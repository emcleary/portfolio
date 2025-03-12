# BUS

A BUS helps a CPU interact with other components of a computer. In a
Game Boy, the BUS is primarily responsible for helping the assembly
instructions read and write data to memory. This does get tricky,
however, since other processing units can affect whether or not data
is accessible to the CPU.

My solution was to use a Mediator design pattern. This helps CPU and
instructions read and write without directly interacting with other
processing units. A simple example from ```BUS::write``` in the
[code](bus.cpp) is the following
```C++
} else if ((0x8000 <= addr) && (addr < 0xA000)) {
    if (!IO::lcds.is_mode_xfer())
        CORE::vram.write(addr, value);
    return;
```
The inner if check is to see if the PPU is currently reading that
memory (pixel data) to display on the screen. It is not safe for the
CPU to edit that data during that time, hence it skips writing
anything.

Overall the BUS might look simple, but it was very important for the
development process. These checks for transferring data were easy to
implement once found and understood. Reading and writing to and from
the APU (audio processing) channels is also very complex, but those
discussions are omitted here.
