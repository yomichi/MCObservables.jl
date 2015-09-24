now_str() = "[" * LibC.strftime("%T", time()) * "] "
inform(str...) = info(now_str(), map(string, str)...)
