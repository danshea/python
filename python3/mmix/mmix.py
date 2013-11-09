#!/usr/bin/env python3

'''
Author: djs
Date: 2012-07-06
Description: mmix emulator - attempts to implement the mmix architecture as
             defined by Donald Knuth.
             
             **************************************************
             *                 WORK IN PROGRESS               *
             **************************************************
'''

class MMIX(object):
    # opcodes
    opcodes = { 'TRAP':0x00,'FCMP':0x01,'FUN':0x02,'FEQL':0x03,'FADD':0x04,'FIX':0x05,'FSUB':0x06,'FIXU':0x07,
                'FLOT':0x08,'FLOTI':0x09,'FLOTU':0x0a,'FLOTUI':0x0b,'SFLOT':0x0c,'SFLOTI':0x0d,'SFLOTU':0x0e,'SFLOTUI':0x0f,
                'FMUL':0x10,'FCMPE':0x11,'FUNE':0x12,'FEQLE':0x13,'FDIV':0x14,'FSQRT':0x15,'FREM':0x16,'FINT':0x17,
                'MUL':0x18,'MULI':0x19,'MULU':0x1a,'MULUI':0x1b,'DIV':0x1c,'DIVI':0x1d,'DIVU':0x1e,'DIVUI':0x1f,
                'ADD':0x20,'ADDI':0x21,'ADDU':0x22,'ADDUI':0x23,'SUB':0x24,'SUBI':0x25,'SUBU':0x26,'SUBUI':0x27,
                '2ADDU':0x28,'2ADDUI':0x29,'4ADDU':0x2a,'4ADDUI':0x2b,'8ADDU':0x2c,'8ADDUI':0x2d,'16ADDU':0x2e,'16ADDUI':0x2f,
                'CMP':0x30,'CMPI':0x31,'CMPU':0x32,'CMPUI':0x33,'NEG':0x34,'NEGI':0x35,'NEGU':0x36,'NEGUI':0x37,
                'SL':0x38,'SLI':0x39,'SLU':0x3a,'SLUI':0x3b,'SR':0x3c,'SRI':0x3d,'SRU':0x3e,'SRUI':0x3f,
                'BN':0x40,'BNB':0x41,'BZ':0x42,'BZB':0x43,'BP':0x44,'BPB':0x45,'BOD':0x46,'BODB':0x47,
                'BNN':0x48,'BNNB':0x49,'BNZ':0x4a,'BNZB':0x4b,'BNP':0x4c,'BNPB':0x4d,'BEV':0x4e,'BEVB':0x4f,
                'PBN':0x50,'PBNB':0x51,'PBZ':0x52,'PBZB':0x53,'PBP':0x54,'PBPB':0x55,'PBOD':0x56,'PBODB':0x57,
                'PBNN':0x58,'PBNNB':0x59,'PBNZ':0x5a,'PBNZB':0x5b,'PBNP':0x5c,'PBNPB':0x5d,'PBEV':0x5e,'PBEVB':0x5f,
                'CSN':0x60,'CSNI':0x61,'CSZ':0x62,'CSZI':0x63,'CSP':0x64,'CSPI':0x65,'CSOD':0x66,'CSODI':0x67,
                'CSNN':0x68,'CSNNI':0x69,'CSNZ':0x6a,'CSNZI':0x6b,'CSNP':0x6c,'CSNPI':0x6d,'CSEV':0x6e,'CSEVI':0x6f,
                'ZSN':0x70,'ZSNI':0x71,'ZSZ':0x72,'ZSZI':0x73,'ZSP':0x74,'ZSPI':0x75,'ZSOD':0x76,'ZSODI':0x77,
                'ZSNN':0x78,'ZSNNI':0x79,'ZSNZ':0x7a,'ZSNZI':0x7b,'ZSNP':0x7c,'ZSNPI':0x7d,'ZSEV':0x7e,'ZSEVI':0x7f,
                'LDB':0x80,'LDBI':0x81,'LDBU':0x82,'LDBUI':0x83,'LDW':0x84,'LDWI':0x85,'LDWU':0x86,'LDWUI':0x87,
                'LDT':0x88,'LDTI':0x89,'LDTU':0x8a,'LDTUI':0x8b,'LDO':0x8c,'LDOI':0x8d,'LDOU':0x8e,'LDOUI':0x8f,
                'LDSF':0x90,'LDSFI':0x91,'LDHT':0x92,'LDHTI':0x93,'CSWAP':0x94,'CSWAPI':0x95,'LDUNC':0x96,'LDUNCI':0x97,
                'LDVTS':0x98,'LDVTSI':0x99,'PRELD':0x9a,'PRELDI':0x9b,'PREGO':0x9c,'PREGOI':0x9d,'GO':0x9e,'GOI':0x9f,
                'STB':0xa0,'STBI':0xa1,'STBU':0xa2,'STBUI':0xa3,'STW':0xa4,'STWI':0xa5,'STWU':0xa6,'STWUI':0xa7,
                'STT':0xa8,'STTI':0xa9,'STTU':0xaa,'STTUI':0xab,'STO':0xac,'STOI':0xad,'STOU':0xae,'STOUI':0xaf,
                'STSF':0xb0,'STSFI':0xb1,'STHT':0xb2,'STHTI':0xb3,'STCO':0xb4,'STCOI':0xb5,'STUNC':0xb6,'STUNCI':0xb7,
                'SYNCD':0xb8,'SYNCDI':0xb9,'PREST':0xba,'PRESTI':0xbb,'SYNCID':0xbc,'SYNCIDI':0xbd,'PUSHGO':0xbe,'PUSHGOI':0xbf,
                'OR':0xc0,'ORI':0xc1,'ORN':0xc2,'ORNI':0xc3,'NOR':0xc4,'NORI':0xc5,'XOR':0xc6,'XORI':0xc7,
                'AND':0xc8,'ANDI':0xc9,'ANDN':0xca,'ANDNI':0xcb,'NAND':0xcc,'NANDI':0xcd,'NXOR':0xce,'NXORI':0xcf,
                'BDIF':0xd0,'BDIFI':0xd1,'WDIF':0xd2,'WDIFI':0xd3,'TDIF':0xd4,'TDIFI':0xd5,'ODIF':0xd6,'ODIFI':0xd7,
                'MUX':0xd8,'MUXI':0xd9,'SADD':0xda,'SADDI':0xdb,'MOR':0xdc,'MORI':0xdd,'MXOR':0xde,'MXORI':0xdf,
                'SETH':0xe0,'SETMH':0xe1,'SETML':0xe2,'SETL':0xe3,'INCH':0xe4,'INCMH':0xe5,'INCML':0xe6,'INCL':0xe7,
                'ORH':0xe8,'ORMH':0xe9,'ORML':0xea,'ORL':0xeb,'ANDNH':0xec,'ANDNMH':0xed,'ANDNML':0xee,'ANDNL':0xef,
                'JMP':0xf0,'JMPB':0xf1,'PUSHJ':0xf2,'PUSHJB':0xf3,'GETA':0xf4,'GETAB':0xf5,'PUT':0xf6,'PUTI':0xf7,
                'POP':0xf8,'RESUME':0xf9,'SAVE':0xfa,'UNSAVE':0xfb,'SYNC':0xfc,'SWYM':0xfd,'GET':0xfe,'TRIP':0xff,
                }

    # define valid sizes
    BYTE  = 8
    WYDE  = 16
    TETRA = 32
    OCTA  = 64

    # define valie modes
    UNSIGNED = 0
    SIGNED   = 1
    FLOAT    = 2

    def __init__(self, memsize=1024):      
        # The 256 General purpose registers
        self.registers = dict(zip(['${0:d}'.format(i) for i in range(256)], ['0'*self.OCTA for i in range(256)]))

        # The 32 Special purpose registers
        self.sregisters = dict(zip(['rA', 'rB', 'rC', 'rD', 'rE', 'rF', 'rG', 'rH', 'rI', 'rJ', 'rK', 'rL', 'rM',
                                    'rN', 'rO', 'rP', 'rQ', 'rR', 'rS', 'rT', 'rU', 'rV', 'rW', 'rX', 'rY', 'rZ',
                                    'rBB', 'rTT', 'rXX', 'rYY', 'rZZ'], ['0'*self.OCTA for i in range(32)]))

        # Memory address space
        self.memsize = memsize
        self.memory = ['0'*self.BYTE for i in range(memsize)]

    '''
    Signed numbers using two's complement notation. (The number's one's complement plus one)
    If the leading bit is a 1, we subtract 2**n to get the integer corresponding to an n-bit number
    in this notation.
    '''

    def sbvalue(self, num, size=BYTE):
        '''
        sbvalue returns the binary value of num
        for BYTE, WYDE, TETRA, OCTA values.
        '''
        offset = 2**(size-1)

        if -offset > num or num > offset-1:
            raise OverflowError('Signed value {0} can not be represented in {1} bits. Valid range is {2} <= n <= {3}'.format(num, size, -offset, offset-1))
        
        if num > 0:
            return bin(num)[2:].zfill(size)
        else:
            return '1'+bin(num+offset)[2:].zfill(size-1)

    def sdvalue(self, num, size=BYTE):
        '''
        sdvalue returns the decimal value of the signed binary
        num of size BYTE, WYDE, TETRA, OCTA
        '''
        if len(num) != size:
            raise ValueError('Argument size mismatch. Received {0} bits, expecting {1}'.format(len(num), size))
        offset = 2**(size-1)
        if num[0] == '1':
            return int(num[1:],2) - offset
        else:
            return int(num[1:],2)
        
    '''
    Unsigned numbers
    '''
    def ubvalue(self, num, size=BYTE):
        '''
        ubvalue returns the binary value of num
        for BYTE, WYDE, TETRA, OCTA values.
        '''
        if num < 0:
            raise ValueError('Unsigned numbers can not be < 0')

        if num > 2**size-1:
            raise OverflowError('Unsigned value {0} can not be represented in {1} bits. Valid range is {2} <= n <= {3}'.format(num, size, 0, 2**size-1))
        
        return bin(num)[2:].zfill(size)

    def udvalue(self, num, size=BYTE):
        '''
        udvalue returns the decimal value of the unsigned binary
        num of size BYTE, WYDE, TETRA, OCTA
        '''
        if len(num) != size:
            raise ValueError('Argument size mismatch. Received {0} bits, expecting {1}'.format(len(num), size))
        
        return int(num,2)
    
    '''
    IEEE 754 Floating Point Representation

    leftmost bit is the sign bit
    next  11 bits for the exponent
    final 52 bits for the mantissa

    1 : 11 : 52 = 64 bit floating point number
    
    1 :  8 : 23 = 32 bit floating point number
    
    Let N be the number we wish to represent.
    exponent e is defined as largest value of 2**e less than N
    mantissa m is defined as 1 - (N / 2**e)
    This is because the one is not included in the mantissa for normal floats
    i.e. 1.m becomes m where m will be kept as an unsigned 23-bit number
    '''

    def fdvalue(self, num, size=TETRA):
        '''
        fdvalue returns the decimal value of the float
        for TETRA, OCTA values of size
        '''
        if size not in [self.TETRA, self.OCTA]:
            raise ValueError('size argument must be of size TETRA ({1}) or OCTA ({1})'.format(TETRA, OCTA))

        if size == self.TETRA:
            esize = 8
        if size == self.OCTA:
            esize = 11

        bias = 2**esize // 2 - 1
        sign = num[0]
        exponent = int(num[1:esize+1],2) - bias
        fraction = num[esize+1:]
        fraction = 1 + sum([int(fraction[i])*2**-(i+1) for i in range(len(fraction))])
        if sign == '1':
            return -(fraction * 2**exponent)
        else:
            return fraction * 2**exponent
        
    def fbvalue(self, num, size=TETRA):
        '''
        fbvalue return a tuple (sign, exponent, fraction) that represents the ieee
        floating point specification.
        '''
        if type(num) is int or type(num) is float:
            num = float(num)
        else:
            raise ValueError('argument must be float or int')

        if size not in [self.TETRA, self.OCTA]:
            raise ValueError('size argument must be of size TETRA ({1}) or OCTA ({1})'.format(TETRA, OCTA))

        if size == self.TETRA:
            esize = 8
            fsize = 23

        if size == self.OCTA:
            esize = 11
            fsize = 52

        # determine the sign bit
        if num < 0:
            sign = '1'
        else:
            sign = '0'

        # sign on the original value is no longer needed
        num = abs(num)
        
        # split the number into its integer and fraction components
        integer, fraction = str(num).split('.')
        integer = int(integer)
        
        # compute the exponent, find the largest value of n such that
        # 2**n <= num
        exponent = esize - 1
        while 2**exponent > num:
            exponent -= 1

        # compute the mantissa (fraction) first for the integer portion
        int_mantissa = integer / 2**exponent

        # determine the offset for the fraction and then compute the
        # mantissa (fraction) for the fraction portion
        # First, take care of leading zero values
        offset = 0
        f = list(fraction)
        while f[0] == '0':
            f.pop(0)
            offset += 1
        # Now determine the offset based on the length of the remaining value
        offset += len(f)
        
        fraction = int(fraction)
        f_mantissa = fraction / 2**exponent * 10**-offset

        mantissa = int_mantissa + f_mantissa - 1

        # convert the mantissa (fraction) to binary
        b_mantissa = ''
        for i in range(1,fsize+1):
            if mantissa >= 2**-i:
                mantissa -= 2**-i
                b_mantissa += '1'
            else:
                b_mantissa += '0'

        # apply the bias to the exponent
        bias = 2**esize // 2 - 1
        exponent = exponent + bias

        return (sign, bin(exponent)[2:].zfill(esize), b_mantissa)

    '''
    Handling memory alignment.
    '''
    
    def align(self, address, size=OCTA):
        '''
        align returns the memory location aligned to byte, wyde, tetra, octa boundaries
        where k = index, we see that:
        M[k]   = M[k]
        M2[2k] = M2[2k+1] = M[2k]M[k+1]
        M4[4k] = ... = M4[4k+3] = M[4k]...M[4k+3]
        M8[8k] = ... = M8[8k+7] = M[8k]...M[8k+7]
        '''
        if address > len(self.memory):
            raise ValueError('Address must be of the range {0} <= n <= {1}'.format(0, len(self.memory)))

        if size not in [self.BYTE, self.WYDE, self.TETRA, self.OCTA]:
            raise ValueError('Alignment size must be of size BYTE({0}) WYDE({1}) TETRA({2}) OCTA({3})'.format(self.BYTE, self.WYDE, self.TETRA, self.OCTA))

        return address - (address % (size//self.BYTE))

    def compute_address(self, Y, Z):
        if Y not in self.registers:
            raise ValueError('{0} is not a valid register.'.format(Y))
        if Z not in self.registers:
            raise ValueError('{0} is not a valid register.'.format(Z))
        address = (self.udvalue(self.registers[Y], self.OCTA) + self.udvalue(self.registers[Z], self.OCTA)) % 2**self.OCTA
        if address >= self.memsize:
            raise ValueError('{0} exceeds maximum memory address size of {1}.'.format(address, self.memsize-1))
        return address
    
    '''
    Loading values.
    '''    
    def load(self, target, address, size=OCTA, mode=UNSIGNED):
        '''
        load into target value from address aligned to size for the given mode
        valid size alignments are BYTE, WYDE, TETRA, OCTA
        valid modes are UNSIGNED, SIGNED, FLOAT
        '''
        if size not in [self.BYTE, self.WYDE, self.TETRA, self.OCTA]:
            raise ValueError('Alignment size must be of size BYTE({0}) WYDE({1}) TETRA({2}) OCTA({3})'.format(self.BYTE, self.WYDE, self.TETRA, self.OCTA))
        if mode not in [self.UNSIGNED, self.SIGNED, self.FLOAT]:
            raise ValueError('Mode must be of UNSIGNED({0}) SIGNED({1}) FLOAT({2})'.format(self.UNSIGNED, self.SIGNED, self.FLOAT))

        start = self.align(address, size)
        end = start + size//self.BYTE
        memval = ''.join(self.memory[start:end])
        if mode == self.UNSIGNED:
            # register is zero filled for unsigned values smaller than an OCTA
            self.registers[target] = memval.zfill(self.OCTA)
        if mode == self.SIGNED:
            # sign bit left fills the register
            self.registers[target] = memval[0]*(self.OCTA-size)+memval
        if mode == self.FLOAT:
            # sign bit left fills the register
            self.registers[target] = memval[0]*(self.OCTA-size)+memval
                             
    '''
    Storing values.
    '''
    def store(self, target, address, size=OCTA, mode=UNSIGNED):
        '''
        store into address, aligned for the given mode, the value in target
        valid size alignments are BYTE, WYDE, TETRA, OCTA
        valid modes are UNSIGNED, SIGNED, FLOAT
        '''
        if size not in [self.BYTE, self.WYDE, self.TETRA, self.OCTA]:
            raise ValueError('Alignment size must be of size BYTE({0}) WYDE({1}) TETRA({2}) OCTA({3})'.format(self.BYTE, self.WYDE, self.TETRA, self.OCTA))
        if mode not in [self.UNSIGNED, self.SIGNED, self.FLOAT]:
            raise ValueError('Mode must be of UNSIGNED({0}) SIGNED({1}) FLOAT({2})'.format(self.UNSIGNED, self.SIGNED, self.FLOAT))

        val = self.registers[target]
        start = self.align(address, size)
        end = start + size//self.BYTE
        val = val[-size:]
        if mode == self.UNSIGNED:
            self.memory[start:end] = [val[i*self.BYTE:i*self.BYTE+self.BYTE] for i in range(0, size//self.BYTE)]
        if mode == self.SIGNED:
            # Behaviour is identical to UNSGINED, but overflow can occur.
            self.memory[start:end] = [val[i*self.BYTE:i*self.BYTE+self.BYTE] for i in range(0, size//self.BYTE)]
        if mode == self.FLOAT:
            self.memory[start:end] = [val[i*self.BYTE:i*self.BYTE+self.BYTE] for i in range(0, size//self.BYTE)]
