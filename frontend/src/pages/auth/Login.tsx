/**
 * Login Page
 */
import React, { useEffect } from 'react'
import { Form, Input, Button, Card, Typography, message, Alert } from 'antd'
import { UserOutlined, LockOutlined } from '@ant-design/icons'
import { useNavigate, Link } from 'react-router-dom'
import { authStore } from '@/store/authStore'
import styles from './Auth.module.css'

const { Title, Text } = Typography

export const Login: React.FC = () => {
  const navigate = useNavigate()
  const { login, isLoading, error, clearError, isAuthenticated } = authStore()
  const [form] = Form.useForm()

  useEffect(() => {
    // Redirect if already authenticated
    if (isAuthenticated) {
      navigate('/dashboard')
    }
  }, [isAuthenticated, navigate])

  useEffect(() => {
    // Clear error when component unmounts
    return () => clearError()
  }, [clearError])

  const handleSubmit = async (values: { username: string; password: string }) => {
    try {
      await login(values)
      message.success('Login successful!')
      navigate('/dashboard')
    } catch (error: any) {
      message.error(error.message || 'Login failed. Please check your credentials.')
    }
  }

  return (
    <div className={styles.authContainer}>
      <Card className={styles.authCard} bordered={false}>
        <div className={styles.authHeader}>
          <Title level={2} className={styles.authTitle}>
            Welcome to NGSmodule
          </Title>
          <Text type="secondary">
            Enterprise Bioinformatics Workstation
          </Text>
        </div>

        {error && (
          <Alert
            message="Login Error"
            description={error}
            type="error"
            closable
            onClose={clearError}
            style={{ marginBottom: 24 }}
          />
        )}

        <Form
          form={form}
          name="login"
          onFinish={handleSubmit}
          autoComplete="off"
          size="large"
          layout="vertical"
        >
          <Form.Item
            name="username"
            rules={[
              { required: true, message: 'Please input your username or email!' },
            ]}
          >
            <Input
              prefix={<UserOutlined />}
              placeholder="Username or Email"
              autoComplete="username"
            />
          </Form.Item>

          <Form.Item
            name="password"
            rules={[{ required: true, message: 'Please input your password!' }]}
          >
            <Input.Password
              prefix={<LockOutlined />}
              placeholder="Password"
              autoComplete="current-password"
            />
          </Form.Item>

          <Form.Item>
            <Button
              type="primary"
              htmlType="submit"
              loading={isLoading}
              block
              size="large"
            >
              Log in
            </Button>
          </Form.Item>
        </Form>

        <div className={styles.authFooter}>
          <Text type="secondary">
            Don't have an account?{' '}
            <Link to="/register" className={styles.authLink}>
              Register now
            </Link>
          </Text>
        </div>
      </Card>
    </div>
  )
}

export default Login
