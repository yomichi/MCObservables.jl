
now_str() = "[" * strftime("%T", time()) * "] "

inform(str...) = info(now_str(), map(string, str)...)
