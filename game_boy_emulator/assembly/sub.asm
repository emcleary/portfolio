INCLUDE "hardware.inc"

SECTION "Header", ROM0[$100]
	jp EntryPoint
	ds $150 - @, 0 ; Make room for the header

EntryPoint:
	ld a, $80
	sub a, a
	ld a, $08
	sub a, a

	ld a, $00
	ld b, $81
	sub a, b

	ld a, $00
	ld b, $80
	ld hl, $8080
	ld [hl], b
	sub a, $80
	sub a, [hl]
	stop
