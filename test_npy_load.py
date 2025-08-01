import os, sys, argparse, logging, random, time
import numpy as np

if __name__ == '__main__':
    
    logging.basicConfig(
        # force=True,
        level=logging.INFO,
        stream=sys.stdout,
        format='### INFO - %(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logging.info(f"{sys.argv[0]} starting up")

    kmc_file = "out.npy"
    # enc_file = "encoded_kmc.txt"

    kmc_bit_obj = np.load(kmc_file)
    kmc_bit_obj_shape = kmc_bit_obj.shape
    # kmc_bit_obj_2d = np.array(kmc_bit_obj.tolist())
    # kmc_bit_obj_2d = kmc_bit_obj_2d[kmc_bit_obj_2d[:, 0].argsort()]
    # kmc_bit_obj = kmc_bit_obj_2d
    print(kmc_bit_obj_shape)
    print(kmc_bit_obj[:10])
    
    first = kmc_bit_obj[0][0]
    second = kmc_bit_obj[9][0]

    b = 16*2
    all_lo = (4**(b//2)-1) // 3

    st_xor = np.bitwise_xor(first, second)
    hamming_dist = np.bitwise_count(np.bitwise_and(np.bitwise_or(np.bitwise_right_shift(st_xor, 1), st_xor), all_lo))
    print(hamming_dist)
    logging.info("Script done.")


            
