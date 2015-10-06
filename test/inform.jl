now_str() = "[" * string(now()) * "] "
inform(str...) = info(now_str(), map(string, str)...)
