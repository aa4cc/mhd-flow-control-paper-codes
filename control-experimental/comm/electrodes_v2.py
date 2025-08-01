import usb.backend.libusb1
import usb.core
import usb.backend
import usb.util

class ElectrodeUSB(object):
    """
    MSG is constructed as:
    [MSG_TYPE, DATA_LEN, DATA]
    -> DATA or DATA_LEN are not mandatory if not needed
    -> only MSG_TYPE and DATA can be send, when the DATA_LEN is always same
    """
    VENDOR_ID = 0x0483
    PRODUCT_ID = 0x5740
    TIMEOUT = 100
    BULK_IN_ADDR = 0x81
    BULK_OUT_ADDR = 0x1
    # --------------------- MSG TYPES -----------------
    SET_LED = 0
    SET_CHANNELS = 1
    SET_ALL_CHANNELS = 2
    CLEAR_CHANNELS = 3
    MEASURE_CURRENT = 4
    MEASURE_ALL_CURRENTS = 5
    GET_FIRMWARE_VERSION = 6
    # --------------- RESPONSE MSG TYPES -----------
    MEASURED_DATA = 7
    FIRMWARE_VERSION = 8
    
    def __init__(self):
        self.el = None
        try:
            backend = usb.backend.libusb1.get_backend(find_library=lambda x: "/opt/homebrew/lib/libusb-1.0.0.dylib")
            self.el = usb.core.find(backend=backend, idVendor=self.VENDOR_ID, idProduct = self.PRODUCT_ID)
        except ValueError as e:
            print(e)
            raise Exception(f"""{e} \n \n[[If 'Device not found' and you recently uploaded code to the board, try plugging and unplogging the UBS
                                 If 'No backend availible' check 'debug.txt']]""")
        if self.el is None:
            raise Exception("[[Device not found]]")
        cfg = usb.util.find_descriptor(self.el, bConfigurationValue=1)
        self.el.set_configuration(cfg)
    
    def read(self, data_len:int, as_str:bool = False):
        # returns accepted data as [uint8] list or as string
        res = self.el.read(self.BULK_IN_ADDR, data_len, self.TIMEOUT)
        if as_str and res is not None:
            return "".join(chr(i) for i in res)
        return res

    def write(self, data:bytes):
        # returns number of bytes written
        return self.el.write(self.BULK_OUT_ADDR, data, self.TIMEOUT)
    
    def set_led(self, value):
        # set leds to either: value = 0 - OFF, 1 - ON, 255 - toggle
        if(value in [0, 1, 255]):
            msg = [self.SET_LED, value]
            self.write(bytes(msg))
        else:
            print("[[|SET LED| - Invalid LED value]]")
        
    def toggle_led(self):
        self.set_led(255)
        
    def set_channels(self, channel_indexes: '[int]', channel_pwms: '[int]'):
        """
        Sets pwm values in 'channel_pwms' to electrodes specified in 'channel_indexes'
        pwm in channel_pwms[0] corresponds to channel_indexes[0]
        pwms can be from 0-255 -> 0 = OFF, 255 = 100% ON
        channel_indexes from 1 to 12
        """
        if (len(channel_pwms) != len(channel_indexes)) or (len(channel_pwms) > 12) or (len(channel_pwms) < 1):
            print("[[|SET CHANNELS| - wrong input data size]]")
            return
        if any((x < 1 or x > 12)for x in channel_indexes):
            print("[[|SET CHANNELS| - wrong channels indexes <1 or >12]]")
            return
        if any((x < 0 or x > 255)for x in channel_pwms):
            print("[[|SET CHANNELS| - wrong channels pwms <0 or >255]]")
            return
        msg = [self.SET_CHANNELS, 2 * len(channel_indexes), *channel_indexes, *channel_pwms]
        self.write(bytes(msg))
    
    def set_all_channels(self, channel_pwms: '[int]'):
        """
        Sets pwm values in 'channel_pwms' to electrodes
        len(channels_pwms) should be 12 and channel_pwms[0] = pwm for electrode 1, pwms[1] = electrode 2 ...
        pwms can be from 0-255 -> 0 = OFF, 255 = 100% ON
        """
        if (len(channel_pwms) != 12):
            print("[[|SET ALL CHANNELS| - wrong input data size]]")
            return
        if any((x < 0 or x > 255)for x in channel_pwms):
            print("[[|SET ALL CHANNELS| - wrong channels pwms <0 or >255]]")
            return

        msg = [self.SET_ALL_CHANNELS, 12, *channel_pwms]
        self.write(bytes(msg))
    
    def clear_channels(self):
        """
        Sets all channels to GND (OFF)
        """
        msg = [self.CLEAR_CHANNELS]
        self.write(bytes(msg))
    
    def measure_current(self, channels_indexes: '[int]'):
        """
        Returns measured current flowing through specified channels_indexes.
        | currently returns only measured voltage in format 0-255 = 0-3,3V to get current - calibration will be needed|
        """
        raise NotImplementedError("Not working yet")
        if (len(channels_indexes) < 1 or len(channels_indexes) > 12):
            print("[[|MEASURE CURRENT| - wrong input data size]]")
            return
        if any((x < 1 or x > 12) for x in channels_indexes):
            print("[[|MEASURE CURRENT| - wrong channels indexes <1 or >12]]")
            return
        msg = [self.MEASURE_CURRENT, len(channels_indexes), *channels_indexes]
        self.write(bytes(msg))
        currents = self.read(len(channels_indexes) + 2)
        if currents[0] != self.MEASURED_DATA:
            print("[[|MEASURE CURRENT| - wrong msg type returned - this should not happen :( ]]")
            print(currents)
            return
        return currents[2:]
    
    def measure_all_currents(self):
        """
        Measures current flowing through all channels - not yet working
        """
        raise NotImplementedError("Not working yet")
        msg = [self.MEASURE_ALL_CURRENTS, 12, *[i for i in range(1, 13)]]
        self.write(bytes(msg))
        currents = self.read(14) # expecting 14 bytes
        if currents[0] != self.MEASURED_DATA:
            print("[[|MEASURE ALL CURRENTS| - wrong msg type returned - this should not happen :( ]]")
            print(currents)
            return
        return currents[2:]

    def get_firmware_version(self):
        # FIRMWARE VERSION IN THIS FORMAT: [x.x.xx] eg [1.0.06] len = 6 bytes (total = 8 bytes)
        msg = [self.GET_FIRMWARE_VERSION]
        self.write(bytes(msg))
        version = self.read(8, True)
        return version[2:]
        

"""
This is without usb protocol - not yet implemented, propably wont be
import serial
import serial.tools.list_ports

BOARD_USB_NAME = "STM32"
USB_BAUDRATE = 4_000_000 # 4Mbit/s - STM32 should support up to 12Mbit/s
class ElectrodeComm(object):
    def __init__(self, port:str=None, baudrate=115_200, use_usb:bool = True):
        # try to find port named BOARD_USB_NAME if port is None and use_usb = True
        self.use_usb = use_usb
        self.port = None
        if port is None:
            port = self.find_usb()
            use_usb = True
            if port is None:
                raise Exception("[[No port provided and Electrode board not connected!]]")
        baud_rate = baudrate if not self.use_usb else USB_BAUDRATE
        try:
            self.port = serial.Serial(port, baud_rate, 8, parity=serial.PARITY_EVEN, stopbits=1)
        except serial.SerialException as e:
            raise Exception(f"Error while openning serial port: {e}")

    def find_usb(self):
        for port in serial.tools.list_ports.comports():
            if port.description == BOARD_USB_NAME:
                return port.device
        return None    
    
    def open(self):
        if not self.port.is_open:
            self.port.open()
        
    def close(self):
        if self.port.is_open:
            self.port.close()
    
    def __del__(self):
        if (self.port is not None and self.port.is_open):
            self.port.close()
    
"""     
"""
if __name__ == "__main__":
    import time
    comm = ElectrodeUSB()
    version = comm.get_firmware_version()
    print(f"[[FIRMWARE VERSION: {version}]]")

    comm.set_all_channels([10]*12)
    time.sleep(2)
    print("Clear_channels")
    comm.clear_channels()
    #print("Measure current on channel 12")
    #currents = comm.measure_current([12])
    #print(f"Measured currents: {currents}")
    #currents = comm.measure_all_currents()
    #print(f"Measured currents: {currents}")
    time.sleep(1)
    
    mapped_states = [12, 10, 7, 5, 2, 4, 1, 3, 6, 8, 11, 9]
    for _ in range(5):
        for i in range(12):
            state = [0] * 12
            state[mapped_states[i] - 1] = 255
            comm.set_all_channels(state)
            time.sleep(0.02)
    comm.clear_channels()
    time.sleep(2)
    state = [0] * 12
    for i in range(12):
        state[mapped_states[i] - 1] = round(((i / 11) * 255))
    print(state)
    comm.set_all_channels(state)
    time.sleep(5)
    comm.clear_channels()
"""