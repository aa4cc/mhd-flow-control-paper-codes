import usb.core
import usb.util
import usb.backend.libusb1
import usb.backend

class CurrentUSB(object):
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
    # --------------- RESPONSE MSG TYPES -----------
    
    # --------------- MAX PWM VALUE
    MAX_PWM = 360
    
    def __init__(self):
        self.el = None
        try:
            backend = usb.backend.libusb1.get_backend(find_library=lambda x: "/opt/homebrew/lib/libusb.dylib")
            self.el = usb.core.find(backend=backend, idVendor=self.VENDOR_ID, idProduct = self.PRODUCT_ID)
        except ValueError as e:
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
        #print(list(data))
        return self.el.write(self.BULK_OUT_ADDR, data, self.TIMEOUT)
    
    def toggle_led(self):
        self.set_led(255)
    
    def set_led(self, value):
        # set leds to either: value = 0 - OFF, 1 - ON, 255 - toggle
        if(value in [0, 1, 255]):
            msg = [self.SET_LED, value]
            self.write(bytes(msg))
        else:
            print("[[|SET LED| - Invalid LED value]]")
    
    def _deconstruct_value(self, value):
        lower = value & 255
        upper = (value >> 8) & 255
        return [upper, lower]
    
    
    def set_channels(self, channel_indexes, channel_currents):
        """
        channel_indexes = 1 to 12
        channel_current = 0.0 to 1.0 (sets current sinking to 0.0 Amp, 1.0 = 1 Amp) or 2.0 HIGH side off 3.0 HIGH side on
        """
        if (len(channel_currents) != len(channel_indexes)) or (len(channel_currents) > 12) or (len(channel_currents) < 1):
            print("[[|SET CHANNELS| - wrong input data size]]")
            return
        if any((x < 1 or x > 12)for x in channel_indexes):
            print("[[|SET CHANNELS| - wrong channels indexes <1 or >12]]")
            return
        if any(((x < 0.0 or x > 1.0) and x != 2.0 and x != 3.0)for x in channel_currents):
            print("[[|SET CHANNELS| - wrong channels currents <0.0 or >1.0]]")
            return
        data = []
        for val in channel_currents:
            value = round(self.MAX_PWM * val)
            data += self._deconstruct_value(value)
        msg = [self.SET_CHANNELS, len(channel_currents), *channel_indexes, *data]
        self.write(bytes(msg))
        
        
"""
if __name__ == "__main__":
    import time
    print("init comm")
    comm = CurrentUSB()
    print("toggle LED")
    comm.set_led(255)
    time.sleep(1)
    print("toggle LED")
    comm.set_led(255)
    time.sleep(1)
    #print("set channel 1 to 0.5A set channel 2 to HIGH")
    #comm.set_channels([1, 2], [0.5, 3.0])
    
    for r in range(10):
        for i in range(1, 13):
            if(r % 2 == 0):
                idx = i
                val = 3.0
            else:
                idx = 13 - i
                val = 2.0
            comm.set_channels([idx], [val])
            time.sleep(1/30)
        time.sleep(0.1)
    time.sleep(1)
    
    
    for r in range(10):
        for i in range(1, 13):
            if(r % 2 == 0):
                idx = i
                val = 1.0
            else:
                idx = 13 - i
                val = 0.0
            comm.set_channels([idx], [val])
            time.sleep(1/30)
        time.sleep(0.1)
    time.sleep(1)
    
    #comm.set_channels([1, 2], [0.0, 2.0])
"""