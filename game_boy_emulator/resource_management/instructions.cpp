#include "instructions.hpp"

#include <memory>

#include "arguments.hpp"
#include "instr.hpp"
#include "../globals.hpp"


std::unique_ptr<Instruction> instr[256];
std::unique_ptr<Instruction> instr_CB[256];


Instruction& get_instruction() {
    uint16_t key = CORE::bus.read();
    IO::timer.add_cycles(1);

    if (key == 0xcb) {
        key = CORE::bus.read();
        IO::timer.add_cycles(1);
        return *instr_CB[key];
    } else {
        return *instr[key];
    }
}


REG::A rA;
REG::B rB;
REG::C rC;
REG::D rD;
REG::E rE;
REG::F rF;
REG::H rH;
REG::L rL;

REG::AF rAF;
REG::BC rBC;
REG::DE rDE;
REG::HL rHL;
REG::SP rSP;
REG::PC rPC;

ADDR::C aC;
ADDR::BC aBC;
ADDR::DE aDE;
ADDR::HL aHL;
ADDR::HLi aHLi;
ADDR::HLd aHLd;

ADDR::N8 aN8;
ADDR::N16 aN16;
ADDR::N16_U16 aN16_U16;

BYTES::E8 bE8;
BYTES::N8 bN8;
BYTES::N16 bN16;

FLAG::Z fZ;
FLAG::NZ fNZ;
FLAG::C fC;
FLAG::NC fNC;


auto fill_vectors = []() {
    instr[0x00] = std::make_unique<NOP>();
    instr[0x01] = std::make_unique<LD_R16_N16>(rBC, bN16);
    instr[0x02] = std::make_unique<LD_AR16_R8>(aBC, rA);
    instr[0x03] = std::make_unique<INC_R16>(rBC);
    instr[0x04] = std::make_unique<INC_R8>(rB);
    instr[0x05] = std::make_unique<DEC_R8>(rB);
    instr[0x06] = std::make_unique<LD_R8_N8>(rB, bN8);
    instr[0x07] = std::make_unique<RLCA>();
    instr[0x08] = std::make_unique<LD_AN16_R16>(aN16_U16, rSP);
    instr[0x09] = std::make_unique<ADD_R16_R16>(rHL, rBC);
    instr[0x0A] = std::make_unique<LD_R8_AR16>(rA, aBC);
    instr[0x0B] = std::make_unique<DEC_R16>(rBC);
    instr[0x0C] = std::make_unique<INC_R8>(rC);
    instr[0x0D] = std::make_unique<DEC_R8>(rC);
    instr[0x0E] = std::make_unique<LD_R8_N8>(rC, bN8);
    instr[0x0F] = std::make_unique<RRCA>();

    instr[0x10] = std::make_unique<STOP>();
    instr[0x11] = std::make_unique<LD_R16_N16>(rDE, bN16);
    instr[0x12] = std::make_unique<LD_AR16_R8>(aDE, rA);
    instr[0x13] = std::make_unique<INC_R16>(rDE);
    instr[0x14] = std::make_unique<INC_R8>(rD);
    instr[0x15] = std::make_unique<DEC_R8>(rD);
    instr[0x16] = std::make_unique<LD_R8_N8>(rD, bN8);
    instr[0x17] = std::make_unique<RLA>();
    instr[0x18] = std::make_unique<JR_E8>(bE8);
    instr[0x19] = std::make_unique<ADD_R16_R16>(rHL, rDE);
    instr[0x1A] = std::make_unique<LD_R8_AR16>(rA, aDE);
    instr[0x1B] = std::make_unique<DEC_R16>(rDE);
    instr[0x1C] = std::make_unique<INC_R8>(rE);
    instr[0x1D] = std::make_unique<DEC_R8>(rE);
    instr[0x1E] = std::make_unique<LD_R8_N8>(rE, bN8);
    instr[0x1F] = std::make_unique<RRA>();

    instr[0x20] = std::make_unique<JR_F_E8>(fNZ, bE8);
    instr[0x21] = std::make_unique<LD_R16_N16>(rHL, bN16);
    instr[0x22] = std::make_unique<LD_AR16_R8>(aHLi, rA);
    instr[0x23] = std::make_unique<INC_R16>(rHL);
    instr[0x24] = std::make_unique<INC_R8>(rH);
    instr[0x25] = std::make_unique<DEC_R8>(rH);
    instr[0x26] = std::make_unique<LD_R8_N8>(rH, bN8);
    instr[0x27] = std::make_unique<DAA>();
    instr[0x28] = std::make_unique<JR_F_E8>(fZ, bE8);
    instr[0x29] = std::make_unique<ADD_R16_R16>(rHL, rHL);
    instr[0x2A] = std::make_unique<LD_R8_AR16>(rA, aHLi);
    instr[0x2B] = std::make_unique<DEC_R16>(rHL);
    instr[0x2C] = std::make_unique<INC_R8>(rL);
    instr[0x2D] = std::make_unique<DEC_R8>(rL);
    instr[0x2E] = std::make_unique<LD_R8_N8>(rL, bN8);
    instr[0x2F] = std::make_unique<CPL>();

    instr[0x30] = std::make_unique<JR_F_E8>(fNC, bE8);
    instr[0x31] = std::make_unique<LD_R16_N16>(rSP, bN16);
    instr[0x32] = std::make_unique<LD_AR16_R8>(aHLd, rA);
    instr[0x33] = std::make_unique<INC_R16>(rSP);
    instr[0x34] = std::make_unique<INC_AR16>(aHL);
    instr[0x35] = std::make_unique<DEC_AR16>(aHL);
    instr[0x36] = std::make_unique<LD_AR16_N8>(aHL, bN8);
    instr[0x37] = std::make_unique<SCF>();
    instr[0x38] = std::make_unique<JR_F_E8>(fC, bE8);
    instr[0x39] = std::make_unique<ADD_R16_R16>(rHL, rSP);
    instr[0x3A] = std::make_unique<LD_R8_AR16>(rA, aHLd);
    instr[0x3B] = std::make_unique<DEC_R16>(rSP);
    instr[0x3C] = std::make_unique<INC_R8>(rA);
    instr[0x3D] = std::make_unique<DEC_R8>(rA);
    instr[0x3E] = std::make_unique<LD_R8_N8>(rA, bN8);
    instr[0x3F] = std::make_unique<CCF>();

    instr[0x40] = std::make_unique<LD_R8_R8>(rB, rB);
    instr[0x41] = std::make_unique<LD_R8_R8>(rB, rC);
    instr[0x42] = std::make_unique<LD_R8_R8>(rB, rD);
    instr[0x43] = std::make_unique<LD_R8_R8>(rB, rE);
    instr[0x44] = std::make_unique<LD_R8_R8>(rB, rH);
    instr[0x45] = std::make_unique<LD_R8_R8>(rB, rL);
    instr[0x46] = std::make_unique<LD_R8_AR16>(rB, aHL);
    instr[0x47] = std::make_unique<LD_R8_R8>(rB, rA);
    instr[0x48] = std::make_unique<LD_R8_R8>(rC, rB);
    instr[0x49] = std::make_unique<LD_R8_R8>(rC, rC);
    instr[0x4A] = std::make_unique<LD_R8_R8>(rC, rD);
    instr[0x4B] = std::make_unique<LD_R8_R8>(rC, rE);
    instr[0x4C] = std::make_unique<LD_R8_R8>(rC, rH);
    instr[0x4D] = std::make_unique<LD_R8_R8>(rC, rL);
    instr[0x4E] = std::make_unique<LD_R8_AR16>(rC, aHL);
    instr[0x4F] = std::make_unique<LD_R8_R8>(rC, rA);

    instr[0x50] = std::make_unique<LD_R8_R8>(rD, rB);
    instr[0x51] = std::make_unique<LD_R8_R8>(rD, rC);
    instr[0x52] = std::make_unique<LD_R8_R8>(rD, rD);
    instr[0x53] = std::make_unique<LD_R8_R8>(rD, rE);
    instr[0x54] = std::make_unique<LD_R8_R8>(rD, rH);
    instr[0x55] = std::make_unique<LD_R8_R8>(rD, rL);
    instr[0x56] = std::make_unique<LD_R8_AR16>(rD, aHL);
    instr[0x57] = std::make_unique<LD_R8_R8>(rD, rA);
    instr[0x58] = std::make_unique<LD_R8_R8>(rE, rB);
    instr[0x59] = std::make_unique<LD_R8_R8>(rE, rC);
    instr[0x5A] = std::make_unique<LD_R8_R8>(rE, rD);
    instr[0x5B] = std::make_unique<LD_R8_R8>(rE, rE);
    instr[0x5C] = std::make_unique<LD_R8_R8>(rE, rH);
    instr[0x5D] = std::make_unique<LD_R8_R8>(rE, rL);
    instr[0x5E] = std::make_unique<LD_R8_AR16>(rE, aHL);
    instr[0x5F] = std::make_unique<LD_R8_R8>(rE, rA);

    instr[0x60] = std::make_unique<LD_R8_R8>(rH, rB);
    instr[0x61] = std::make_unique<LD_R8_R8>(rH, rC);
    instr[0x62] = std::make_unique<LD_R8_R8>(rH, rD);
    instr[0x63] = std::make_unique<LD_R8_R8>(rH, rE);
    instr[0x64] = std::make_unique<LD_R8_R8>(rH, rH);
    instr[0x65] = std::make_unique<LD_R8_R8>(rH, rL);
    instr[0x66] = std::make_unique<LD_R8_AR16>(rH, aHL);
    instr[0x67] = std::make_unique<LD_R8_R8>(rH, rA);
    instr[0x68] = std::make_unique<LD_R8_R8>(rL, rB);
    instr[0x69] = std::make_unique<LD_R8_R8>(rL, rC);
    instr[0x6A] = std::make_unique<LD_R8_R8>(rL, rD);
    instr[0x6B] = std::make_unique<LD_R8_R8>(rL, rE);
    instr[0x6C] = std::make_unique<LD_R8_R8>(rL, rH);
    instr[0x6D] = std::make_unique<LD_R8_R8>(rL, rL);
    instr[0x6E] = std::make_unique<LD_R8_AR16>(rL, aHL);
    instr[0x6F] = std::make_unique<LD_R8_R8>(rL, rA);

    instr[0x70] = std::make_unique<LD_AR16_R8>(aHL, rB);
    instr[0x71] = std::make_unique<LD_AR16_R8>(aHL, rC);
    instr[0x72] = std::make_unique<LD_AR16_R8>(aHL, rD);
    instr[0x73] = std::make_unique<LD_AR16_R8>(aHL, rE);
    instr[0x74] = std::make_unique<LD_AR16_R8>(aHL, rH);
    instr[0x75] = std::make_unique<LD_AR16_R8>(aHL, rL);
    instr[0x76] = std::make_unique<HALT>();
    instr[0x77] = std::make_unique<LD_AR16_R8>(aHL, rA);
    instr[0x78] = std::make_unique<LD_R8_R8>(rA, rB);
    instr[0x79] = std::make_unique<LD_R8_R8>(rA, rC);
    instr[0x7A] = std::make_unique<LD_R8_R8>(rA, rD);
    instr[0x7B] = std::make_unique<LD_R8_R8>(rA, rE);
    instr[0x7C] = std::make_unique<LD_R8_R8>(rA, rH);
    instr[0x7D] = std::make_unique<LD_R8_R8>(rA, rL);
    instr[0x7E] = std::make_unique<LD_R8_AR16>(rA, aHL);
    instr[0x7F] = std::make_unique<LD_R8_R8>(rA, rA);

    instr[0x80] = std::make_unique<ADD_R8_R8>(rA, rB);
    instr[0x81] = std::make_unique<ADD_R8_R8>(rA, rC);
    instr[0x82] = std::make_unique<ADD_R8_R8>(rA, rD);
    instr[0x83] = std::make_unique<ADD_R8_R8>(rA, rE);
    instr[0x84] = std::make_unique<ADD_R8_R8>(rA, rH);
    instr[0x85] = std::make_unique<ADD_R8_R8>(rA, rL);
    instr[0x86] = std::make_unique<ADD_R8_AR16>(rA, aHL);
    instr[0x87] = std::make_unique<ADD_R8_R8>(rA, rA);
    instr[0x88] = std::make_unique<ADC_R8_R8>(rA, rB);
    instr[0x89] = std::make_unique<ADC_R8_R8>(rA, rC);
    instr[0x8A] = std::make_unique<ADC_R8_R8>(rA, rD);
    instr[0x8B] = std::make_unique<ADC_R8_R8>(rA, rE);
    instr[0x8C] = std::make_unique<ADC_R8_R8>(rA, rH);
    instr[0x8D] = std::make_unique<ADC_R8_R8>(rA, rL);
    instr[0x8E] = std::make_unique<ADC_R8_AR16>(rA, aHL);
    instr[0x8F] = std::make_unique<ADC_R8_R8>(rA, rA);

    instr[0x90] = std::make_unique<SUB_R8_R8>(rA, rB);
    instr[0x91] = std::make_unique<SUB_R8_R8>(rA, rC);
    instr[0x92] = std::make_unique<SUB_R8_R8>(rA, rD);
    instr[0x93] = std::make_unique<SUB_R8_R8>(rA, rE);
    instr[0x94] = std::make_unique<SUB_R8_R8>(rA, rH);
    instr[0x95] = std::make_unique<SUB_R8_R8>(rA, rL);
    instr[0x96] = std::make_unique<SUB_R8_AR16>(rA, aHL);
    instr[0x97] = std::make_unique<SUB_R8_R8>(rA, rA);
    instr[0x98] = std::make_unique<SBC_R8_R8>(rA, rB);
    instr[0x99] = std::make_unique<SBC_R8_R8>(rA, rC);
    instr[0x9A] = std::make_unique<SBC_R8_R8>(rA, rD);
    instr[0x9B] = std::make_unique<SBC_R8_R8>(rA, rE);
    instr[0x9C] = std::make_unique<SBC_R8_R8>(rA, rH);
    instr[0x9D] = std::make_unique<SBC_R8_R8>(rA, rL);
    instr[0x9E] = std::make_unique<SBC_R8_AR16>(rA, aHL);
    instr[0x9F] = std::make_unique<SBC_R8_R8>(rA, rA);

    instr[0xA0] = std::make_unique<AND_R8_R8>(rA, rB);
    instr[0xA1] = std::make_unique<AND_R8_R8>(rA, rC);
    instr[0xA2] = std::make_unique<AND_R8_R8>(rA, rD);
    instr[0xA3] = std::make_unique<AND_R8_R8>(rA, rE);
    instr[0xA4] = std::make_unique<AND_R8_R8>(rA, rH);
    instr[0xA5] = std::make_unique<AND_R8_R8>(rA, rL);
    instr[0xA6] = std::make_unique<AND_R8_AR16>(rA, aHL);
    instr[0xA7] = std::make_unique<AND_R8_R8>(rA, rA);
    instr[0xA8] = std::make_unique<XOR_R8_R8>(rA, rB);
    instr[0xA9] = std::make_unique<XOR_R8_R8>(rA, rC);
    instr[0xAA] = std::make_unique<XOR_R8_R8>(rA, rD);
    instr[0xAB] = std::make_unique<XOR_R8_R8>(rA, rE);
    instr[0xAC] = std::make_unique<XOR_R8_R8>(rA, rH);
    instr[0xAD] = std::make_unique<XOR_R8_R8>(rA, rL);
    instr[0xAE] = std::make_unique<XOR_R8_AR16>(rA, aHL);
    instr[0xAF] = std::make_unique<XOR_R8_R8>(rA, rA);

    instr[0xB0] = std::make_unique<OR_R8_R8>(rA, rB);
    instr[0xB1] = std::make_unique<OR_R8_R8>(rA, rC);
    instr[0xB2] = std::make_unique<OR_R8_R8>(rA, rD);
    instr[0xB3] = std::make_unique<OR_R8_R8>(rA, rE);
    instr[0xB4] = std::make_unique<OR_R8_R8>(rA, rH);
    instr[0xB5] = std::make_unique<OR_R8_R8>(rA, rL);
    instr[0xB6] = std::make_unique<OR_R8_AR16>(rA, aHL);
    instr[0xB7] = std::make_unique<OR_R8_R8>(rA, rA);
    instr[0xB8] = std::make_unique<CP_R8_R8>(rA, rB);
    instr[0xB9] = std::make_unique<CP_R8_R8>(rA, rC);
    instr[0xBA] = std::make_unique<CP_R8_R8>(rA, rD);
    instr[0xBB] = std::make_unique<CP_R8_R8>(rA, rE);
    instr[0xBC] = std::make_unique<CP_R8_R8>(rA, rH);
    instr[0xBD] = std::make_unique<CP_R8_R8>(rA, rL);
    instr[0xBE] = std::make_unique<CP_R8_AR16>(rA, aHL);
    instr[0xBF] = std::make_unique<CP_R8_R8>(rA, rA);

    instr[0xC0] = std::make_unique<RET_F>(fNZ);
    instr[0xC1] = std::make_unique<POP>(rBC);
    instr[0xC2] = std::make_unique<JP_F_N16>(fNZ, bN16);
    instr[0xC3] = std::make_unique<JP_N16>(bN16);
    instr[0xC4] = std::make_unique<CALL_F_N16>(fNZ, bN16);
    instr[0xC5] = std::make_unique<PUSH>(rBC);
    instr[0xC6] = std::make_unique<ADD_R8_N8>(rA, bN8);
    instr[0xC7] = std::make_unique<RST>(0x00);
    instr[0xC8] = std::make_unique<RET_F>(fZ);
    instr[0xC9] = std::make_unique<RET>();
    instr[0xCA] = std::make_unique<JP_F_N16>(fZ, bN16);
    instr[0xCB] = std::make_unique<DUMMY>(0xCB);
    instr[0xCC] = std::make_unique<CALL_F_N16>(fZ, bN16);
    instr[0xCD] = std::make_unique<CALL_N16>(bN16);
    instr[0xCE] = std::make_unique<ADC_R8_N8>(rA, bN8);
    instr[0xCF] = std::make_unique<RST>(0x08);

    instr[0xD0] = std::make_unique<RET_F>(fNC);
    instr[0xD1] = std::make_unique<POP>(rDE);
    instr[0xD2] = std::make_unique<JP_F_N16>(fNC, bN16);
    instr[0xD3] = std::make_unique<DUMMY>(0xD3);
    instr[0xD4] = std::make_unique<CALL_F_N16>(fNC, bN16);
    instr[0xD5] = std::make_unique<PUSH>(rDE);
    instr[0xD6] = std::make_unique<SUB_R8_N8>(rA, bN8);
    instr[0xD7] = std::make_unique<RST>(0x10);
    instr[0xD8] = std::make_unique<RET_F>(fC);
    instr[0xD9] = std::make_unique<RETI>();
    instr[0xDA] = std::make_unique<JP_F_N16>(fC, bN16);
    instr[0xDB] = std::make_unique<DUMMY>(0xDB);
    instr[0xDC] = std::make_unique<CALL_F_N16>(fC, bN16);
    instr[0xDD] = std::make_unique<DUMMY>(0xDD);
    instr[0xDE] = std::make_unique<SBC_R8_N8>(rA, bN8);
    instr[0xDF] = std::make_unique<RST>(0x18);

    instr[0xE0] = std::make_unique<LD_AN8_R8>(aN8, rA);
    instr[0xE1] = std::make_unique<POP>(rHL);
    instr[0xE2] = std::make_unique<LD_AR8_R8>(aC, rA);
    instr[0xE3] = std::make_unique<DUMMY>(0xE3);
    instr[0xE4] = std::make_unique<DUMMY>(0xE4);
    instr[0xE5] = std::make_unique<PUSH>(rHL);
    instr[0xE6] = std::make_unique<AND_R8_N8>(rA, bN8);
    instr[0xE7] = std::make_unique<RST>(0x20);
    instr[0xE8] = std::make_unique<ADD_R16_E8>(rSP, bE8);
    instr[0xE9] = std::make_unique<JP_R16>(rHL);
    instr[0xEA] = std::make_unique<LD_AN16_R8>(aN16, rA);
    instr[0xEB] = std::make_unique<DUMMY>(0xEB);
    instr[0xEC] = std::make_unique<DUMMY>(0xEC);
    instr[0xED] = std::make_unique<DUMMY>(0xED);
    instr[0xEE] = std::make_unique<XOR_R8_N8>(rA, bN8);
    instr[0xEF] = std::make_unique<RST>(0x28);

    instr[0xF0] = std::make_unique<LD_R8_AN8>(rA, aN8);
    instr[0xF1] = std::make_unique<POP>(rAF);
    instr[0xF2] = std::make_unique<LD_R8_AR8>(rA, aC);
    instr[0xF3] = std::make_unique<DI>();
    instr[0xF4] = std::make_unique<DUMMY>(0xF4);
    instr[0xF5] = std::make_unique<PUSH>(rAF);
    instr[0xF6] = std::make_unique<OR_R8_N8>(rA, bN8);
    instr[0xF7] = std::make_unique<RST>(0x30);
    instr[0xF8] = std::make_unique<LD_HL_SP_E8>(rHL, rSP, bE8);
    instr[0xF9] = std::make_unique<LD_R16_R16>(rSP, rHL);
    instr[0xFA] = std::make_unique<LD_R8_AN16>(rA, aN16);
    instr[0xFB] = std::make_unique<EI>();
    instr[0xFC] = std::make_unique<DUMMY>(0xFC);
    instr[0xFD] = std::make_unique<DUMMY>(0xFD);
    instr[0xFE] = std::make_unique<CP_R8_N8>(rA, bN8);
    instr[0xFF] = std::make_unique<RST>(0x38);


    instr_CB[0x00] = std::make_unique<RLC_R8>(rB);
    instr_CB[0x01] = std::make_unique<RLC_R8>(rC);
    instr_CB[0x02] = std::make_unique<RLC_R8>(rD);
    instr_CB[0x03] = std::make_unique<RLC_R8>(rE);
    instr_CB[0x04] = std::make_unique<RLC_R8>(rH);
    instr_CB[0x05] = std::make_unique<RLC_R8>(rL);
    instr_CB[0x06] = std::make_unique<RLC_AR16>(aHL);
    instr_CB[0x07] = std::make_unique<RLC_R8>(rA);

    instr_CB[0x08] = std::make_unique<RRC_R8>(rB);
    instr_CB[0x09] = std::make_unique<RRC_R8>(rC);
    instr_CB[0x0A] = std::make_unique<RRC_R8>(rD);
    instr_CB[0x0B] = std::make_unique<RRC_R8>(rE);
    instr_CB[0x0C] = std::make_unique<RRC_R8>(rH);
    instr_CB[0x0D] = std::make_unique<RRC_R8>(rL);
    instr_CB[0x0E] = std::make_unique<RRC_AR16>(aHL);
    instr_CB[0x0F] = std::make_unique<RRC_R8>(rA);

    instr_CB[0x10] = std::make_unique<RL_R8>(rB);
    instr_CB[0x11] = std::make_unique<RL_R8>(rC);
    instr_CB[0x12] = std::make_unique<RL_R8>(rD);
    instr_CB[0x13] = std::make_unique<RL_R8>(rE);
    instr_CB[0x14] = std::make_unique<RL_R8>(rH);
    instr_CB[0x15] = std::make_unique<RL_R8>(rL);
    instr_CB[0x16] = std::make_unique<RL_AR16>(aHL);
    instr_CB[0x17] = std::make_unique<RL_R8>(rA);

    instr_CB[0x18] = std::make_unique<RR_R8>(rB);
    instr_CB[0x19] = std::make_unique<RR_R8>(rC);
    instr_CB[0x1A] = std::make_unique<RR_R8>(rD);
    instr_CB[0x1B] = std::make_unique<RR_R8>(rE);
    instr_CB[0x1C] = std::make_unique<RR_R8>(rH);
    instr_CB[0x1D] = std::make_unique<RR_R8>(rL);
    instr_CB[0x1E] = std::make_unique<RR_AR16>(aHL);
    instr_CB[0x1F] = std::make_unique<RR_R8>(rA);

    instr_CB[0x20] = std::make_unique<SLA_R8>(rB);
    instr_CB[0x21] = std::make_unique<SLA_R8>(rC);
    instr_CB[0x22] = std::make_unique<SLA_R8>(rD);
    instr_CB[0x23] = std::make_unique<SLA_R8>(rE);
    instr_CB[0x24] = std::make_unique<SLA_R8>(rH);
    instr_CB[0x25] = std::make_unique<SLA_R8>(rL);
    instr_CB[0x26] = std::make_unique<SLA_AR16>(aHL);
    instr_CB[0x27] = std::make_unique<SLA_R8>(rA);

    instr_CB[0x28] = std::make_unique<SRA_R8>(rB);
    instr_CB[0x29] = std::make_unique<SRA_R8>(rC);
    instr_CB[0x2A] = std::make_unique<SRA_R8>(rD);
    instr_CB[0x2B] = std::make_unique<SRA_R8>(rE);
    instr_CB[0x2C] = std::make_unique<SRA_R8>(rH);
    instr_CB[0x2D] = std::make_unique<SRA_R8>(rL);
    instr_CB[0x2E] = std::make_unique<SRA_AR16>(aHL);
    instr_CB[0x2F] = std::make_unique<SRA_R8>(rA);

    instr_CB[0x30] = std::make_unique<SWAP_R8>(rB);
    instr_CB[0x31] = std::make_unique<SWAP_R8>(rC);
    instr_CB[0x32] = std::make_unique<SWAP_R8>(rD);
    instr_CB[0x33] = std::make_unique<SWAP_R8>(rE);
    instr_CB[0x34] = std::make_unique<SWAP_R8>(rH);
    instr_CB[0x35] = std::make_unique<SWAP_R8>(rL);
    instr_CB[0x36] = std::make_unique<SWAP_AR16>(aHL);
    instr_CB[0x37] = std::make_unique<SWAP_R8>(rA);

    instr_CB[0x38] = std::make_unique<SRL_R8>(rB);
    instr_CB[0x39] = std::make_unique<SRL_R8>(rC);
    instr_CB[0x3A] = std::make_unique<SRL_R8>(rD);
    instr_CB[0x3B] = std::make_unique<SRL_R8>(rE);
    instr_CB[0x3C] = std::make_unique<SRL_R8>(rH);
    instr_CB[0x3D] = std::make_unique<SRL_R8>(rL);
    instr_CB[0x3E] = std::make_unique<SRL_AR16>(aHL);
    instr_CB[0x3F] = std::make_unique<SRL_R8>(rA);

    instr_CB[0x40] = std::make_unique<BIT_R8>(0, rB);
    instr_CB[0x41] = std::make_unique<BIT_R8>(0, rC);
    instr_CB[0x42] = std::make_unique<BIT_R8>(0, rD);
    instr_CB[0x43] = std::make_unique<BIT_R8>(0, rE);
    instr_CB[0x44] = std::make_unique<BIT_R8>(0, rH);
    instr_CB[0x45] = std::make_unique<BIT_R8>(0, rL);
    instr_CB[0x46] = std::make_unique<BIT_AR16>(0, aHL);
    instr_CB[0x47] = std::make_unique<BIT_R8>(0, rA);

    instr_CB[0x48] = std::make_unique<BIT_R8>(1, rB);
    instr_CB[0x49] = std::make_unique<BIT_R8>(1, rC);
    instr_CB[0x4A] = std::make_unique<BIT_R8>(1, rD);
    instr_CB[0x4B] = std::make_unique<BIT_R8>(1, rE);
    instr_CB[0x4C] = std::make_unique<BIT_R8>(1, rH);
    instr_CB[0x4D] = std::make_unique<BIT_R8>(1, rL);
    instr_CB[0x4E] = std::make_unique<BIT_AR16>(1, aHL);
    instr_CB[0x4F] = std::make_unique<BIT_R8>(1, rA);

    instr_CB[0x50] = std::make_unique<BIT_R8>(2, rB);
    instr_CB[0x51] = std::make_unique<BIT_R8>(2, rC);
    instr_CB[0x52] = std::make_unique<BIT_R8>(2, rD);
    instr_CB[0x53] = std::make_unique<BIT_R8>(2, rE);
    instr_CB[0x54] = std::make_unique<BIT_R8>(2, rH);
    instr_CB[0x55] = std::make_unique<BIT_R8>(2, rL);
    instr_CB[0x56] = std::make_unique<BIT_AR16>(2, aHL);
    instr_CB[0x57] = std::make_unique<BIT_R8>(2, rA);

    instr_CB[0x58] = std::make_unique<BIT_R8>(3, rB);
    instr_CB[0x59] = std::make_unique<BIT_R8>(3, rC);
    instr_CB[0x5A] = std::make_unique<BIT_R8>(3, rD);
    instr_CB[0x5B] = std::make_unique<BIT_R8>(3, rE);
    instr_CB[0x5C] = std::make_unique<BIT_R8>(3, rH);
    instr_CB[0x5D] = std::make_unique<BIT_R8>(3, rL);
    instr_CB[0x5E] = std::make_unique<BIT_AR16>(3, aHL);
    instr_CB[0x5F] = std::make_unique<BIT_R8>(3, rA);

    instr_CB[0x60] = std::make_unique<BIT_R8>(4, rB);
    instr_CB[0x61] = std::make_unique<BIT_R8>(4, rC);
    instr_CB[0x62] = std::make_unique<BIT_R8>(4, rD);
    instr_CB[0x63] = std::make_unique<BIT_R8>(4, rE);
    instr_CB[0x64] = std::make_unique<BIT_R8>(4, rH);
    instr_CB[0x65] = std::make_unique<BIT_R8>(4, rL);
    instr_CB[0x66] = std::make_unique<BIT_AR16>(4, aHL);
    instr_CB[0x67] = std::make_unique<BIT_R8>(4, rA);

    instr_CB[0x68] = std::make_unique<BIT_R8>(5, rB);
    instr_CB[0x69] = std::make_unique<BIT_R8>(5, rC);
    instr_CB[0x6A] = std::make_unique<BIT_R8>(5, rD);
    instr_CB[0x6B] = std::make_unique<BIT_R8>(5, rE);
    instr_CB[0x6C] = std::make_unique<BIT_R8>(5, rH);
    instr_CB[0x6D] = std::make_unique<BIT_R8>(5, rL);
    instr_CB[0x6E] = std::make_unique<BIT_AR16>(5, aHL);
    instr_CB[0x6F] = std::make_unique<BIT_R8>(5, rA);

    instr_CB[0x70] = std::make_unique<BIT_R8>(6, rB);
    instr_CB[0x71] = std::make_unique<BIT_R8>(6, rC);
    instr_CB[0x72] = std::make_unique<BIT_R8>(6, rD);
    instr_CB[0x73] = std::make_unique<BIT_R8>(6, rE);
    instr_CB[0x74] = std::make_unique<BIT_R8>(6, rH);
    instr_CB[0x75] = std::make_unique<BIT_R8>(6, rL);
    instr_CB[0x76] = std::make_unique<BIT_AR16>(6, aHL);
    instr_CB[0x77] = std::make_unique<BIT_R8>(6, rA);

    instr_CB[0x78] = std::make_unique<BIT_R8>(7, rB);
    instr_CB[0x79] = std::make_unique<BIT_R8>(7, rC);
    instr_CB[0x7A] = std::make_unique<BIT_R8>(7, rD);
    instr_CB[0x7B] = std::make_unique<BIT_R8>(7, rE);
    instr_CB[0x7C] = std::make_unique<BIT_R8>(7, rH);
    instr_CB[0x7D] = std::make_unique<BIT_R8>(7, rL);
    instr_CB[0x7E] = std::make_unique<BIT_AR16>(7, aHL);
    instr_CB[0x7F] = std::make_unique<BIT_R8>(7, rA);

    instr_CB[0x80] = std::make_unique<RES_R8>(0, rB);
    instr_CB[0x81] = std::make_unique<RES_R8>(0, rC);
    instr_CB[0x82] = std::make_unique<RES_R8>(0, rD);
    instr_CB[0x83] = std::make_unique<RES_R8>(0, rE);
    instr_CB[0x84] = std::make_unique<RES_R8>(0, rH);
    instr_CB[0x85] = std::make_unique<RES_R8>(0, rL);
    instr_CB[0x86] = std::make_unique<RES_AR16>(0, aHL);
    instr_CB[0x87] = std::make_unique<RES_R8>(0, rA);

    instr_CB[0x88] = std::make_unique<RES_R8>(1, rB);
    instr_CB[0x89] = std::make_unique<RES_R8>(1, rC);
    instr_CB[0x8A] = std::make_unique<RES_R8>(1, rD);
    instr_CB[0x8B] = std::make_unique<RES_R8>(1, rE);
    instr_CB[0x8C] = std::make_unique<RES_R8>(1, rH);
    instr_CB[0x8D] = std::make_unique<RES_R8>(1, rL);
    instr_CB[0x8E] = std::make_unique<RES_AR16>(1, aHL);
    instr_CB[0x8F] = std::make_unique<RES_R8>(1, rA);

    instr_CB[0x90] = std::make_unique<RES_R8>(2, rB);
    instr_CB[0x91] = std::make_unique<RES_R8>(2, rC);
    instr_CB[0x92] = std::make_unique<RES_R8>(2, rD);
    instr_CB[0x93] = std::make_unique<RES_R8>(2, rE);
    instr_CB[0x94] = std::make_unique<RES_R8>(2, rH);
    instr_CB[0x95] = std::make_unique<RES_R8>(2, rL);
    instr_CB[0x96] = std::make_unique<RES_AR16>(2, aHL);
    instr_CB[0x97] = std::make_unique<RES_R8>(2, rA);

    instr_CB[0x98] = std::make_unique<RES_R8>(3, rB);
    instr_CB[0x99] = std::make_unique<RES_R8>(3, rC);
    instr_CB[0x9A] = std::make_unique<RES_R8>(3, rD);
    instr_CB[0x9B] = std::make_unique<RES_R8>(3, rE);
    instr_CB[0x9C] = std::make_unique<RES_R8>(3, rH);
    instr_CB[0x9D] = std::make_unique<RES_R8>(3, rL);
    instr_CB[0x9E] = std::make_unique<RES_AR16>(3, aHL);
    instr_CB[0x9F] = std::make_unique<RES_R8>(3, rA);

    instr_CB[0xA0] = std::make_unique<RES_R8>(4, rB);
    instr_CB[0xA1] = std::make_unique<RES_R8>(4, rC);
    instr_CB[0xA2] = std::make_unique<RES_R8>(4, rD);
    instr_CB[0xA3] = std::make_unique<RES_R8>(4, rE);
    instr_CB[0xA4] = std::make_unique<RES_R8>(4, rH);
    instr_CB[0xA5] = std::make_unique<RES_R8>(4, rL);
    instr_CB[0xA6] = std::make_unique<RES_AR16>(4, aHL);
    instr_CB[0xA7] = std::make_unique<RES_R8>(4, rA);

    instr_CB[0xA8] = std::make_unique<RES_R8>(5, rB);
    instr_CB[0xA9] = std::make_unique<RES_R8>(5, rC);
    instr_CB[0xAA] = std::make_unique<RES_R8>(5, rD);
    instr_CB[0xAB] = std::make_unique<RES_R8>(5, rE);
    instr_CB[0xAC] = std::make_unique<RES_R8>(5, rH);
    instr_CB[0xAD] = std::make_unique<RES_R8>(5, rL);
    instr_CB[0xAE] = std::make_unique<RES_AR16>(5, aHL);
    instr_CB[0xAF] = std::make_unique<RES_R8>(5, rA);

    instr_CB[0xB0] = std::make_unique<RES_R8>(6, rB);
    instr_CB[0xB1] = std::make_unique<RES_R8>(6, rC);
    instr_CB[0xB2] = std::make_unique<RES_R8>(6, rD);
    instr_CB[0xB3] = std::make_unique<RES_R8>(6, rE);
    instr_CB[0xB4] = std::make_unique<RES_R8>(6, rH);
    instr_CB[0xB5] = std::make_unique<RES_R8>(6, rL);
    instr_CB[0xB6] = std::make_unique<RES_AR16>(6, aHL);
    instr_CB[0xB7] = std::make_unique<RES_R8>(6, rA);

    instr_CB[0xB8] = std::make_unique<RES_R8>(7, rB);
    instr_CB[0xB9] = std::make_unique<RES_R8>(7, rC);
    instr_CB[0xBA] = std::make_unique<RES_R8>(7, rD);
    instr_CB[0xBB] = std::make_unique<RES_R8>(7, rE);
    instr_CB[0xBC] = std::make_unique<RES_R8>(7, rH);
    instr_CB[0xBD] = std::make_unique<RES_R8>(7, rL);
    instr_CB[0xBE] = std::make_unique<RES_AR16>(7, aHL);
    instr_CB[0xBF] = std::make_unique<RES_R8>(7, rA);

    instr_CB[0xC0] = std::make_unique<SET_R8>(0, rB);
    instr_CB[0xC1] = std::make_unique<SET_R8>(0, rC);
    instr_CB[0xC2] = std::make_unique<SET_R8>(0, rD);
    instr_CB[0xC3] = std::make_unique<SET_R8>(0, rE);
    instr_CB[0xC4] = std::make_unique<SET_R8>(0, rH);
    instr_CB[0xC5] = std::make_unique<SET_R8>(0, rL);
    instr_CB[0xC6] = std::make_unique<SET_AR16>(0, aHL);
    instr_CB[0xC7] = std::make_unique<SET_R8>(0, rA);

    instr_CB[0xC8] = std::make_unique<SET_R8>(1, rB);
    instr_CB[0xC9] = std::make_unique<SET_R8>(1, rC);
    instr_CB[0xCA] = std::make_unique<SET_R8>(1, rD);
    instr_CB[0xCB] = std::make_unique<SET_R8>(1, rE);
    instr_CB[0xCC] = std::make_unique<SET_R8>(1, rH);
    instr_CB[0xCD] = std::make_unique<SET_R8>(1, rL);
    instr_CB[0xCE] = std::make_unique<SET_AR16>(1, aHL);
    instr_CB[0xCF] = std::make_unique<SET_R8>(1, rA);

    instr_CB[0xD0] = std::make_unique<SET_R8>(2, rB);
    instr_CB[0xD1] = std::make_unique<SET_R8>(2, rC);
    instr_CB[0xD2] = std::make_unique<SET_R8>(2, rD);
    instr_CB[0xD3] = std::make_unique<SET_R8>(2, rE);
    instr_CB[0xD4] = std::make_unique<SET_R8>(2, rH);
    instr_CB[0xD5] = std::make_unique<SET_R8>(2, rL);
    instr_CB[0xD6] = std::make_unique<SET_AR16>(2, aHL);
    instr_CB[0xD7] = std::make_unique<SET_R8>(2, rA);

    instr_CB[0xD8] = std::make_unique<SET_R8>(3, rB);
    instr_CB[0xD9] = std::make_unique<SET_R8>(3, rC);
    instr_CB[0xDA] = std::make_unique<SET_R8>(3, rD);
    instr_CB[0xDB] = std::make_unique<SET_R8>(3, rE);
    instr_CB[0xDC] = std::make_unique<SET_R8>(3, rH);
    instr_CB[0xDD] = std::make_unique<SET_R8>(3, rL);
    instr_CB[0xDE] = std::make_unique<SET_AR16>(3, aHL);
    instr_CB[0xDF] = std::make_unique<SET_R8>(3, rA);

    instr_CB[0xE0] = std::make_unique<SET_R8>(4, rB);
    instr_CB[0xE1] = std::make_unique<SET_R8>(4, rC);
    instr_CB[0xE2] = std::make_unique<SET_R8>(4, rD);
    instr_CB[0xE3] = std::make_unique<SET_R8>(4, rE);
    instr_CB[0xE4] = std::make_unique<SET_R8>(4, rH);
    instr_CB[0xE5] = std::make_unique<SET_R8>(4, rL);
    instr_CB[0xE6] = std::make_unique<SET_AR16>(4, aHL);
    instr_CB[0xE7] = std::make_unique<SET_R8>(4, rA);

    instr_CB[0xE8] = std::make_unique<SET_R8>(5, rB);
    instr_CB[0xE9] = std::make_unique<SET_R8>(5, rC);
    instr_CB[0xEA] = std::make_unique<SET_R8>(5, rD);
    instr_CB[0xEB] = std::make_unique<SET_R8>(5, rE);
    instr_CB[0xEC] = std::make_unique<SET_R8>(5, rH);
    instr_CB[0xED] = std::make_unique<SET_R8>(5, rL);
    instr_CB[0xEE] = std::make_unique<SET_AR16>(5, aHL);
    instr_CB[0xEF] = std::make_unique<SET_R8>(5, rA);

    instr_CB[0xF0] = std::make_unique<SET_R8>(6, rB);
    instr_CB[0xF1] = std::make_unique<SET_R8>(6, rC);
    instr_CB[0xF2] = std::make_unique<SET_R8>(6, rD);
    instr_CB[0xF3] = std::make_unique<SET_R8>(6, rE);
    instr_CB[0xF4] = std::make_unique<SET_R8>(6, rH);
    instr_CB[0xF5] = std::make_unique<SET_R8>(6, rL);
    instr_CB[0xF6] = std::make_unique<SET_AR16>(6, aHL);
    instr_CB[0xF7] = std::make_unique<SET_R8>(6, rA);

    instr_CB[0xF8] = std::make_unique<SET_R8>(7, rB);
    instr_CB[0xF9] = std::make_unique<SET_R8>(7, rC);
    instr_CB[0xFA] = std::make_unique<SET_R8>(7, rD);
    instr_CB[0xFB] = std::make_unique<SET_R8>(7, rE);
    instr_CB[0xFC] = std::make_unique<SET_R8>(7, rH);
    instr_CB[0xFD] = std::make_unique<SET_R8>(7, rL);
    instr_CB[0xFE] = std::make_unique<SET_AR16>(7, aHL);
    instr_CB[0xFF] = std::make_unique<SET_R8>(7, rA);

    return true;
}();
