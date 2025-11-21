/**
 * Register Page
 */
import React, { useEffect } from 'react'
import { Form, Input, Button, Card, Typography, message, Alert } from 'antd'
import { UserOutlined, LockOutlined, MailOutlined, TeamOutlined, IdcardOutlined } from '@ant-design/icons'
import { useNavigate, Link } from 'react-router-dom'
import { authStore } from '@/store/authStore'
import styles from './Auth.module.css'

const { Title, Text } = Typography

export const Register: React.FC = () => {
  const navigate = useNavigate()
  const { register, isLoading, error, clearError, isAuthenticated } = authStore()
  const [form] = Form.useForm()

  useEffect(() => {
    // Redirect if already authenticated
    if (isAuthenticated) {
      navigate('/dashboard')
    }
  }, [isAuthenticated, navigate])

  useEffect(() => {
    return () => clearError()
  }, [clearError])

  const handleSubmit = async (values: any) => {
    try {
      await register(values)
      message.success('Registration successful! Welcome to NGSmodule.')
      navigate('/dashboard')
    } catch (error: any) {
      message.error(error.message || 'Registration failed. Please try again.')
    }
  }

  return (
    <div className={styles.authContainer}>
      <Card className={styles.authCard} bordered={false}>
        <div className={styles.authHeader}>
          <Title level={2} className={styles.authTitle}>
            Create Account
          </Title>
          <Text type="secondary">
            Join NGSmodule Bioinformatics Platform
          </Text>
        </div>

        {error && (
          <Alert
            message="Registration Error"
            description={error}
            type="error"
            closable
            onClose={clearError}
            style={{ marginBottom: 24 }}
          />
        )}

        <Form
          form={form}
          name="register"
          onFinish={handleSubmit}
          autoComplete="off"
          size="large"
          layout="vertical"
        >
          <Form.Item
            name="username"
            rules={[
              { required: true, message: 'Please input your username!' },
              { min: 3, message: 'Username must be at least 3 characters!' },
              { max: 50, message: 'Username must be at most 50 characters!' },
              {
                pattern: /^[a-zA-Z0-9_-]+$/,
                message: 'Username can only contain letters, numbers, underscore and hyphen!'
              },
            ]}
          >
            <Input
              prefix={<UserOutlined />}
              placeholder="Username"
              autoComplete="username"
            />
          </Form.Item>

          <Form.Item
            name="email"
            rules={[
              { required: true, message: 'Please input your email!' },
              { type: 'email', message: 'Please input a valid email!' },
            ]}
          >
            <Input
              prefix={<MailOutlined />}
              placeholder="Email"
              autoComplete="email"
            />
          </Form.Item>

          <Form.Item
            name="password"
            rules={[
              { required: true, message: 'Please input your password!' },
              { min: 8, message: 'Password must be at least 8 characters!' },
            ]}
          >
            <Input.Password
              prefix={<LockOutlined />}
              placeholder="Password (min 8 characters)"
              autoComplete="new-password"
            />
          </Form.Item>

          <Form.Item
            name="confirm_password"
            dependencies={['password']}
            rules={[
              { required: true, message: 'Please confirm your password!' },
              ({ getFieldValue }) => ({
                validator(_, value) {
                  if (!value || getFieldValue('password') === value) {
                    return Promise.resolve()
                  }
                  return Promise.reject(new Error('Passwords do not match!'))
                },
              }),
            ]}
          >
            <Input.Password
              prefix={<LockOutlined />}
              placeholder="Confirm Password"
              autoComplete="new-password"
            />
          </Form.Item>

          <Form.Item name="full_name">
            <Input
              prefix={<IdcardOutlined />}
              placeholder="Full Name (optional)"
              autoComplete="name"
            />
          </Form.Item>

          <Form.Item name="organization">
            <Input
              prefix={<TeamOutlined />}
              placeholder="Organization (optional)"
              autoComplete="organization"
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
              Create Account
            </Button>
          </Form.Item>
        </Form>

        <div className={styles.authFooter}>
          <Text type="secondary">
            Already have an account?{' '}
            <Link to="/login" className={styles.authLink}>
              Log in
            </Link>
          </Text>
        </div>
      </Card>
    </div>
  )
}

export default Register
