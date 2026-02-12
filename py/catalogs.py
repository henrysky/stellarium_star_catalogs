import pathlib
import struct
import pandas as pd
import numpy as np
import tqdm


def decode_star_hip(encoded):
    # Ensure the encoded value is exactly 3 bytes
    if len(encoded) != 3:
        raise ValueError("Encoded value must be 3 bytes.")
    
    # Unpack the 3 bytes into a 24-bit integer (little-endian format)
    combined_value = struct.unpack("<I", encoded.ljust(4, b'\0'))[0]
    
    # Extract the hip (17-bit) and letter_value (5-bit)
    hip = combined_value >> 5
    letter_value = combined_value & 0x1F  # Mask to get the lower 5 bits
    
    return hip, letter_value


def decode_star_hip_bigendian(encoded):
    """Decode hip and componentid from 3 bytes in big endian format"""
    # Ensure the encoded value is exactly 3 bytes
    if len(encoded) != 3:
        raise ValueError("Encoded value must be 3 bytes.")
    
    # Unpack the 3 bytes into a 24-bit integer (big-endian format)
    # Prepend with zero byte to make it 4 bytes for unpacking
    combined_value = struct.unpack(">I", b'\0' + encoded)[0]
    
    # Extract componentid (5-bit, stored in upper bits) and hip (19-bit, stored in lower bits)
    componentid = (combined_value >> 19) & 0x1F  # Get upper 5 bits
    hip = combined_value & 0x7FFFF  # Get lower 19 bits
    
    return hip, componentid


def encode_star_hip_big_endian(hip, componentid):
    """Encode hip and componentid into 3 bytes in big endian format"""
    # Combine hip (17-bit) and componentid (5-bit) into a 24-bit value
    # Sacrifice the top 5 bits for componentid, store hip in the lower 19 bits
    combined_value = ((componentid & 0x1F) << 19) | (hip & 0x7FFFF)
    
    # Pack as big endian 32-bit integer and take the last 3 bytes
    packed = struct.pack(">I", combined_value)
    return packed[1:]  # Return last 3 bytes


def encode_uint32_to_3bytes_big_endian(value):
    """Encode a uint32 value into 3 bytes in big endian format"""
    packed = struct.pack(">I", value)
    return packed[1:]  # Return last 3 bytes


def read_to_dataframe(file: pathlib.Path) -> pd.DataFrame:
    f = open(file, "rb")

    # Read the header (28 bytes)
    header = struct.unpack("6if", f.read(28))
    header = dict(
        zip(
            (
                "Magic",
                "Data Type",
                "Major Version",
                "Minor Version",
                "Level",
                "Magnitude Minimum",
                "Epoch",
            ),
            header,
        )
    )
    n_zones = 20 * 4 ** header["Level"] + 1

    # Unpack the zone information
    zone_info = struct.unpack(f"{n_zones}I", f.read(n_zones * 4))
    max_records = sum(zone_info)

    # Define the header for the star data
    if header["Data Type"] == 0:
        star_header = (
            "source_id",
            "x0",
            "x1",
            "x2",
            "dx0",
            "dx1",
            "dx2",
            "b_v",
            "mag",
            "plx",
            "plx_err",
            "rv",
            "sp_int",
            "otype",
            "hip",
            "componentid",
        )
    elif header["Data Type"] == 1:
        star_header = (
            "source_id",
            "x0",
            "x1",
            "dx0",
            "dx1",
            "b_v",
            "mag",
            "plx",
            "plx_err",
        )
    else:
        star_header = (
            "source_id",
            "x0",
            "x1",
            "b_v",
            "mag",
        )

    # Create a dataframe to store the stars
    df = pd.DataFrame(columns=star_header, index=range(max_records))

    # put header in dataframe for reference with df.attrs
    for i in header:
        df.attrs[i] = header[i]

    for i in tqdm.tqdm(range(max_records), desc="Reading stars"):
        if header["Data Type"] == 0:
            star = struct.unpack("qiiiiiihhHHhHB", f.read(45))
            hip, componentid = decode_star_hip(f.read(3))
            df.loc[i] = star + (hip, componentid, )
        elif header["Data Type"] == 1:
            star = struct.unpack("qiiiihhHH", f.read(32))
            # add to dataframe
            df.loc[i] = star
        else:
            source_id = struct.unpack("q", f.read(8))[0]
            # recover uint32 from 3 bytes
            x0 = struct.unpack("<I", b"\0" + f.read(3))[0]
            x1 = struct.unpack("<I", b"\0" + f.read(3))[0]
            b_v = struct.unpack("b", f.read(1))[0]
            mag = struct.unpack("b", f.read(1))[0]
            df.loc[i] = (source_id, x0, x1, b_v, mag)

    # create zone data array
    zone_data = np.zeros(max_records, dtype=int)
    i = 0
    for idx, n in enumerate(zone_info):
        zone_data[i:i+n] = idx
        i += n
    df["zone"] = zone_data

    return df


def convert_to_big_endian(input_file: pathlib.Path, output_file: pathlib.Path) -> None:
    """
    Convert a binary star catalog file from little endian to big endian format.
    
    Args:
        input_file: Path to the input little endian binary file
        output_file: Path to the output big endian binary file
    """
    with open(input_file, "rb") as f_in, open(output_file, "wb") as f_out:
        # Read and convert the header (28 bytes)
        header_data = f_in.read(28)
        header = struct.unpack("<6if", header_data)  # Read as little endian
        
        # Write header in big endian format
        f_out.write(struct.pack(">6if", *header))
        
        # Extract header information for processing
        header_dict = dict(
            zip(
                (
                    "Magic",
                    "Data Type",
                    "Major Version",
                    "Minor Version",
                    "Level",
                    "Magnitude Minimum",
                    "Epoch",
                ),
                header,
            )
        )
        
        n_zones = 20 * 4 ** header_dict["Level"] + 1
        
        # Read and convert zone information
        zone_data = f_in.read(n_zones * 4)
        zone_info = struct.unpack(f"<{n_zones}I", zone_data)  # Read as little endian
        f_out.write(struct.pack(f">{n_zones}I", *zone_info))  # Write as big endian
        
        max_records = sum(zone_info)
        
        # Convert star data based on data type
        for i in tqdm.tqdm(range(max_records), desc="Converting stars to big endian"):
            if header_dict["Data Type"] == 0:
                # Read 45 bytes of star data + 3 bytes of hip data
                star_data = f_in.read(45)
                hip_data = f_in.read(3)
                
                # Unpack star data from little endian
                star = struct.unpack("<qiiiiiihhHHhHB", star_data)
                
                # Write star data in big endian
                f_out.write(struct.pack(">qiiiiiihhHHhHB", *star))
                
                # Decode hip data and re-encode in big endian
                hip, componentid = decode_star_hip(hip_data)
                f_out.write(encode_star_hip_big_endian(hip, componentid))
                
            elif header_dict["Data Type"] == 1:
                # Read 32 bytes of star data
                star_data = f_in.read(32)
                star = struct.unpack("<qiiiihhHH", star_data)
                
                # Write in big endian
                f_out.write(struct.pack(">qiiiihhHH", *star))
                
            else:
                # Data Type 2: Read individual components
                source_id_data = f_in.read(8)
                x0_data = f_in.read(3)
                x1_data = f_in.read(3)
                b_v_data = f_in.read(1)
                mag_data = f_in.read(1)
                
                # Unpack from little endian
                source_id = struct.unpack("<q", source_id_data)[0]
                x0 = struct.unpack("<I", b"\0" + x0_data)[0]
                x1 = struct.unpack("<I", b"\0" + x1_data)[0]
                b_v = struct.unpack("b", b_v_data)[0]
                mag = struct.unpack("b", mag_data)[0]
                
                # Write in big endian
                f_out.write(struct.pack(">q", source_id))
                f_out.write(encode_uint32_to_3bytes_big_endian(x0))
                f_out.write(encode_uint32_to_3bytes_big_endian(x1))
                f_out.write(struct.pack("b", b_v))
                f_out.write(struct.pack("b", mag))
    
    print(f"Successfully converted {input_file} to big endian format: {output_file}")
